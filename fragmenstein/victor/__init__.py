########################################################################################################################

__doc__ = \
    """
Victor (after Dr Victor Frankenstein) is a class that uses both Fragmenstein (makes blended compounds) and Egor (energy minimises).
This master reanimator keeps a ``.journal`` (logging, class attribute).
And can be called via the class method ``.laboratory`` where he can process multiple compounds at once.

    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################

from ..egor import Egor, pyrosetta
from ..core import Fragmenstein
from ..m_rmsd import mRSMD

from typing import List, Tuple, Dict, Union, Optional, Callable

from rdkit_to_params import Params, Constraints

import os, logging, sys, re, pymol2, warnings, json, unicodedata, requests, shutil

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdmolfiles, rdFMCS, Draw


class Victor:
    """
    * ``smiles`` SMILES string (inputted)
    * ``long_name`` name for files
    * ``ligand_resn`` the residue name for the ligand.
    * ``ligand_resi`` the residue index (PDB) for the ligand.
    * ``covalent_resi`` the residue index (PDB) for the covalent attachment
    * ``covalent_resn`` the residue name for the covalent attachment. For now can only be 'CYS'
    * ``params`` Params instance
    * ``constraint`` Constraint or None depending on if covalent.
    * ``mol`` the molecule
    * ``covalent_definitions`` class attr. that stores for each possible attachment residue (CYS) defs for constraints.
    * ``warhead_definitions`` class attr. that stores warheader info
    * ``journal`` class attr. logging
    * ``work_path`` class attr. where to save stuff

    ``warhead_definitions`` and ``covalent_definitions`` are class attributes that can be modified beforehand to
    allow a new attachment. ``covalent_definitions`` is a list of dictionaries of 'residue', 'smiles', 'names',
    which are needed for the constraint file making. Namely smiles is two atoms and the connection and names is the
    names of each. Cysteine is ``{'residue': 'CYS', 'smiles': '*SC', 'names': ['CONN3', 'SG', 'CB']}``.
    While ``warhead_definitions`` is a list of 'name' (name of warhead for humans),
    'covalent' (the smiles of the warhead, where the zeroth atom is the one attached to the rest),
    'noncovalent' (the warhead unreacted),
    'covalent_atomnames' and 'noncovalent_atomnames' (list of atom names).
    The need for atomnames is actually not for the code but to allow lazy tweaks and analysis downstream
    (say typing in pymol: `show sphere, name CX`).
    Adding a 'constraint' to an entry will apply that constraint.

    """
    hits_path = 'hits'
    work_path = 'output'
    journal = logging.getLogger('Fragmenstein')
    journal.setLevel(logging.DEBUG)

    covalent_definitions = [{'residue': 'CYS', 'smiles': '*SC', 'atomnames': ['CONN3', 'SG', 'CB']}]
    warhead_definitions = [{'name': 'nitrile',
                            'covalent': 'C(=N)*',  # zeroth atom is attached to the rest
                            'covalent_atomnames': ['CX', 'NX', 'CONN1'],
                            'noncovalent': 'C(#N)',  # zeroth atom is attached to the rest
                            'noncovalent_atomnames': ['CX', 'NX']
                            },
                           {'name': 'acrylamide',
                            'covalent': 'C(=O)CC*',  # the N may be secondary etc. so best not do mad substitutions.
                            'covalent_atomnames': ['CZ', 'OZ', 'CY', 'CX', 'CONN1'],
                            # OZ needs to tautomerise & h-bond happily.
                            'noncovalent': 'C(=O)C=C',
                            'noncovalent_atomnames': ['CZ', 'OZ', 'CY', 'CX']},
                           {'name': 'chloroacetamide',
                            'covalent': 'C(=O)C*',  # the N may be secondary etc. so best not do mad substitutions.
                            'covalent_atomnames': ['CY', 'OY', 'CX', 'CONN1'],
                            # OY needs to tautomerise & h-bond happily.
                            'noncovalent': 'C(=O)C[Cl]',
                            'noncovalent_atomnames': ['CY', 'OY', 'CX', 'CLX']
                            },
                           {'name': 'vinylsulfonamide',
                            'covalent': 'S(=O)(=O)CC*',  # the N may be secondary etc. so best not do mad substitutions.
                            'covalent_atomnames': ['SZ', 'OZ1', 'OZ2', 'CY', 'CX', 'CONN1'],  # OZ tauto
                            'noncovalent': 'S(=O)(=O)C=C',
                            'noncovalent_atomnames': ['SZ', 'OZ1', 'OZ2', 'CY', 'CX']
                            }
                           ]
    _connected_names = ('CONN', 'LOWE', 'UPPE', 'CONN1', 'CONN2', 'CONN3', 'LOWER', 'UPPER')

    def __init__(self,
                 smiles: str,
                 hits: List[Chem.Mol],
                 pdb_filename: str,
                 long_name: str = 'ligand',
                 ligand_resn: str = 'LIG',
                 ligand_resi: Union[int, str] = '1B',
                 covalent_resn: str = 'CYS',  # no other option is accepted.
                 covalent_resi: Optional[Union[int, str]] = None,
                 extra_constraint: Union[str] = None,
                 pose_fx: Optional[Callable] = None,
                 ):
        """
        :param smiles: smiles of followup, optionally covalent (_e.g._ ``*CC(=O)CCC``)
        :param hits: list of rdkit molecules
        :param pdb_filename: file of apo structure
        :param long_name: gets used for filenames so will get slugified
        :param ligand_resn: 3 letter code or your choice
        :param ligand_resi: Rosetta-style pose(int) or pdb(str)
        :param covalent_resn: only CYS accepted. if smiles has no * it is ignored
        :param covalent_resi: Rosetta-style pose(int) or pdb(str)
        :param extra_constraint: multiline string of constraints..
        :param pose_fx: a function to call with pose to tweak or change something before minimising.
        """
        # ***** STORE *******
        # entry attributes
        self.long_name = self.slugify(long_name)
        self.smiles = smiles
        self.apo_pdbblock = open(pdb_filename).read()
        self.hits = hits
        self.ligand_resn = ligand_resn.upper()
        self.ligand_resi = ligand_resi
        self.covalent_resn = covalent_resn.upper()
        self.covalent_resi = covalent_resi
        self.extra_constraint = extra_constraint
        self.pose_fx = pose_fx
        # warnings
        with warnings.catch_warnings(record=True) as self._warned:
            try:
                # check they are okay
                if '*' in self.smiles and (self.covalent_resi is None or self.covalent_resn is None):
                    raise ValueError(f'{self.long_name} - is covalent but without known covalent residues')
                    # TODO '*' in self.smiles is bad. user might start with a mol file.
                elif '*' in self.smiles:
                    self.is_covalent = True
                else:
                    self.is_covalent = False
                self._assert_inputs()
                # ***** PARAMS & CONSTRAINT *******
                self.journal.info(f'{self.long_name} - Starting work')
                self._log_warnings()
                # making folder.
                self._make_output_folder()
                # make params
                self.journal.debug(f'{self.long_name} - Starting parameterisation')
                self.params = Params.from_smiles(self.smiles, name=ligand_resn, generic=False)
                self.journal.warning(f'{self.long_name} - CHI HAS BEEN DISABLED')
                self.params.CHI.data = []  # TODO fix chi
                self.mol = self.params.mol
                self._log_warnings()
                # deal with covalent and non covalent separately
                if self.is_covalent:
                    self.journal.debug(f'{self.long_name} - is covalent.')
                    self.constraint = self._fix_covalent()
                    if extra_constraint:
                        self.constraint.custom_constraint += self.extra_constraint
                    attachment = self._get_attachment_from_pdbblock()
                else:
                    self.journal.debug(f'{self.long_name} - is not covalent.')
                    if extra_constraint:
                        self.constraint = Constraints.mock()
                        self.constraint.custom_constraint += self.extra_constraint
                    else:
                        self.constraint = None
                    attachment = None
                self._log_warnings()
                # ***** FRAGMENSTEIN *******
                # make fragmenstein
                self.journal.debug(f'{self.long_name} - Starting fragmenstein')
                self.fragmenstein = Fragmenstein(self.mol, self.hits, attachment=attachment)
                self.unminimised_pdbblock = self._place_fragmenstein()
                self.constraint.custom_constraint += self._make_coordinate_constraints()
                self._checkpoint_bravo()
                # save stuff
                params_file, holo_file, constraint_file = self._checkpoint_alpha()
                # ***** EGOR *******
                self.journal.debug(f'{self.long_name} - setting up Egor')
                self.egor = Egor.from_pdbblock(pdbblock=self.unminimised_pdbblock,
                                               params_file=params_file,
                                               constraint_file=constraint_file,
                                               ligand_residue=self.ligand_resi,
                                               key_residues=[self.covalent_resi])
                # user custom code.
                if self.pose_fx is not None:
                    self.journal.debug(f'{self.long_name} - running custom pose mod.')
                    self.pose_fx(self.egor.pose)
                # storing a roundtrip
                self.unminimised_pdbblock = self.egor.pose2str()
                # minimise
                self.journal.debug(f'{self.long_name} - Egor minimising')
                self.egor.minimise()
                self._checkpoint_charlie()
                self.journal.debug(f'{self.long_name} - Completed')
            except Exception as err:
                self.journal.exception(f'{self.long_name} — {err.__class__.__name__}: {err}')

    # =================== Init called methods ==========================================================================

    def slugify(self, name: str):
        return re.sub(r'[\W_.-]+', '-', name)

    def _make_output_folder(self):
        path = os.path.join(self.work_path, self.long_name)
        if not os.path.exists(self.work_path):
            os.mkdir(self.work_path)
        if not os.path.exists(path):
            os.mkdir(path)
        else:
            self.journal.warning(f'{self.long_name} - Folder {path} exists.')

    def _assert_inputs(self):
        assert len(self.ligand_resn) == 3, f'{self.long_name} - {self.ligand_resn} is not 3 char long.'
        assert len(self.hits), f'{self.long_name} - No hits to use to construct'
        assert self.ligand_resn != 'UNL', f'{self.long_name} - It cannot be UNL as it s the unspecified resn in rdkit'
        if self.covalent_resn and len(
                [d for d in self.covalent_definitions if d['residue'] == self.covalent_resn]) == 0:
            raise ValueError(f'{self.long_name} - Unrecognised type {self.covalent_resn}')

    def _fix_covalent(self):
        self.journal.debug(f'{self.long_name} - fixing for covalent')
        # to make life easier for analysis, CX is the attachment atom, CY is the one before it.
        for war_def in self.warhead_definitions:
            warhead = Chem.MolFromSmiles(war_def['covalent'])
            if self.params.mol.HasSubstructMatch(warhead):
                self.params.rename_by_template(warhead, war_def['covalent_atomnames'])
                cov_def = [d for d in self.covalent_definitions if d['residue'] == self.covalent_resn][0]
                self.journal.debug(f'{self.long_name} - has a {war_def["name"]}')
                cons = Constraints(smiles=(war_def['covalent'], cov_def['smiles']),
                                   names=[*war_def['covalent_atomnames'], *cov_def['atomnames']],
                                   ligand_res=self.ligand_resi,
                                   target_res=self.covalent_resi)
                # user added constraint
                if 'constraint' in war_def:
                    cons.custom_constraint = war_def['constraint']
                return cons
        else:
            raise ValueError(f'{self.long_name} - Unsure what the warhead is.')

    def _make_coordinate_constraints(self):
        lines = []
        origins = self.fragmenstein.origin_from_mol(self.fragmenstein.positioned_mol)
        std = self.fragmenstein.stdev_from_mol(self.fragmenstein.positioned_mol)
        conf = self.fragmenstein.positioned_mol.GetConformer()
        for i in range(self.fragmenstein.positioned_mol.GetNumAtoms()):
            if origins[i]:
                atom = self.fragmenstein.positioned_mol.GetAtomWithIdx(i)
                if atom.GetSymbol() in self._connected_names:
                    continue
                pos = conf.GetAtomPosition(i)
                lines.append(f'CoordinateConstraint {atom.GetPDBResidueInfo().GetName()} {self.ligand_resi} '+ \
                             f'CA {self.covalent_resi} '+ \
                             f'{pos.x} {pos.y} {pos.z} HARMONIC 0 {std[i] + 1}\n')
        return ''.join(lines)

    def _get_attachment_from_pdbblock(self) -> Chem.Mol:
        """
        Yes, yes, I see the madness in using pymol to get an atom for rdkit to make a pose for pyrosetta.
        """
        self.journal.debug(f'{self.long_name} - getting attachemnt atom')
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(self.apo_pdbblock, 'prot')
            name = self.constraint.target_con_name.strip()
            resi = re.match('(\d+)', str(self.constraint.target_res)).group(1)
            try:
                chain = re.match('\D', str(self.constraint.target_res)).group(1)
                pdb = pymol.cmd.get_pdbstr(f'resi {resi} and name {name} and chain {chain}')
            except:
                pdb = pymol.cmd.get_pdbstr(f'resi {resi} and name {name}')
            return Chem.MolFromPDBBlock(pdb)

    def _place_fragmenstein(self):
        l_resi, l_chain = re.match('(\d+)(\D?)', str(self.ligand_resi)).groups()
        p_resi, p_chain = re.match('(\d+)(\D?)', str(self.covalent_resi)).groups()
        if not p_chain:
            p_chain = 'A'
        if not l_chain:
            l_chain = 'B'
        self.journal.debug(f'{self.long_name} - placing fragmenstein')
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(self.apo_pdbblock, 'apo')
            # distort positions
            pos_mol = Chem.MolToPDBBlock(self.fragmenstein.positioned_mol)
            pymol.cmd.read_pdbstr(pos_mol, 'scaffold')
            pymol.cmd.alter('scaffold', f'resi="{l_resi}"')
            pymol.cmd.alter('scaffold', f'chain="{l_chain}"')
            pymol.cmd.remove('name R')  # no dummy atoms!
            for c in self._connected_names:
                pymol.cmd.remove(f'name {c}')  # no conns
            pymol.cmd.remove('resn UNL')  # no unmatched stuff.
            pdbblock = pymol.cmd.get_pdbstr('*')
            pymol.cmd.delete('*')
        if self.is_covalent:
            cx = self.params.pad_name(self.params.CONNECT[0].atom_name)
            return f'LINK         SG  {self.covalent_resn} {p_chain} {p_resi: >3}                {cx} {self.ligand_resn} {l_chain} {l_resi: >3}     1555   1555  1.8\n' + pdbblock
        else:
            return pdbblock

    @classmethod
    def copy_names(cls, acceptor_mol, donor_mol):
        mcs = rdFMCS.FindMCS([acceptor_mol, donor_mol],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrder,
                             ringMatchesRingOnly=True)
        common = Chem.MolFromSmarts(mcs.smartsString)
        pos_match = acceptor_mol.positioned_mol.GetSubstructMatch(common)
        pdb_match = donor_mol.GetSubstructMatch(common)
        for m, p in zip(pos_match, pdb_match):
            ma = acceptor_mol.GetAtomWithIdx(m)
            pa = donor_mol.GetAtomWithIdx(p)
            assert ma.GetSymbol() == pa.GetSymbol(), 'The indices do not align! ' + \
                                                     f'{ma.GetIdx()}:{ma.GetSymbol()} vs. ' + \
                                                     f'{pa.GetIdx()}:{pa.GetSymbol()}'
            ma.SetMonomerInfo(pa.GetPDBResidueInfo())

    def _log_warnings(self):
        if len(self._warned):
            for w in self._warned:
                self.journal.warning(f'{self.long_name} - {w.message} ({w.category})')
            self._warned.clear()

    # =================== Logging ======================================================================================

    def _checkpoint_alpha(self):
        self._log_warnings()
        #  saving params
        self.journal.debug(f'{self.long_name} - saving params')
        params_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params')
        self.params.dump(params_file)
        # saving holo
        self.journal.debug(f'{self.long_name} - saving holo (unmimised)')
        holo_file = os.path.join(self.work_path, self.long_name, self.long_name + '.holo_unminimised.pdb')
        with open(holo_file, 'w') as w:
            w.write(self.unminimised_pdbblock)
        # saving constraint
        if self.constraint:
            self.journal.debug(f'{self.long_name} - saving constraint')
            constraint_file = os.path.join(self.work_path, self.long_name, self.long_name + '.con')
            self.constraint.dump(constraint_file)
        else:
            constraint_file = ''
        # saving hits (without copying)
        for h, hit in enumerate(self.hits):
            if hit.HasProp("_Name") and hit.GetProp("_Name").strip():
                name = hit.GetProp("_Name")
            else:
                name = f'hit{h}'
            hfile = os.path.join(self.work_path, self.long_name, f'{name}.pdb')
            Chem.MolToPDBFile(hit, hfile)
            mfile = os.path.join(self.work_path, self.long_name, f'{name}.mol')
            Chem.MolToMolFile(hit, mfile, kekulize=False)
        # saving params template
        params_template_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params_template.pdb')
        Chem.MolToPDBFile(self.params.mol, params_template_file)
        params_template_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params_template.mol')
        Chem.MolToMolFile(self.params.mol, params_template_file)
        # checking all is in order
        self.journal.debug(f'{self.long_name} - checking params file works')
        ptest_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params_test.pdb')
        Params.params_to_pose(params_file, self.params.NAME).dump_pdb(ptest_file)
        self._log_warnings()
        return params_file, holo_file, constraint_file

    def _checkpoint_bravo(self):
        self._log_warnings()
        self.journal.debug(f'{self.long_name} - saving mols from fragmenstein')
        scaffold_file = os.path.join(self.work_path, self.long_name, self.long_name + '.scaffold.mol')
        Chem.MolToMolFile(self.fragmenstein.scaffold, scaffold_file, kekulize=False)
        chimera_file = os.path.join(self.work_path, self.long_name, self.long_name + '.chimera.mol')
        Chem.MolToMolFile(self.fragmenstein.chimera, chimera_file, kekulize=False)
        pos_file = os.path.join(self.work_path, self.long_name, self.long_name + '.positioned.mol')
        Chem.MolToMolFile(self.fragmenstein.positioned_mol, pos_file, kekulize=False)
        frag_file = os.path.join(self.work_path, self.long_name, self.long_name + '.fragmenstein.json')
        with open(frag_file, 'w') as w:
            json.dump({'smiles': self.smiles,
                       'origin': self.fragmenstein.origin_from_mol(self.fragmenstein.positioned_mol),
                       'stdev': self.fragmenstein.stdev_from_mol(self.fragmenstein.positioned_mol)},
                      w)
        self._log_warnings()
        # unminimised_pdbblock will be saved by egor (round trip via pose)

    def _checkpoint_charlie(self):
        self._log_warnings()
        self.journal.debug(f'{self.long_name} - saving pose from egor')
        min_file = os.path.join(self.work_path, self.long_name, self.long_name + '.holo_minimised.pdb')
        self.egor.pose.dump_pdb(min_file)
        self.journal.debug(f'{self.long_name} - calculating Gibbs')
        energy = self.egor.ligand_score()
        # recover bonds
        self.journal.debug(f'{self.long_name} - making ligand only')
        lig_file = os.path.join(self.work_path, self.long_name, self.long_name + '.minimised.mol')
        ligand = self.egor.mol_from_pose()
        template = AllChem.DeleteSubstructs(self.params.mol, Chem.MolFromSmiles('*'))
        ligand = AllChem.AssignBondOrdersFromTemplate(template, ligand)
        self.journal.debug(f'{self.long_name} - calculating mRMSD')
        mrsmd = mRSMD.from_other_annotated_mols(ligand, self.hits, self.fragmenstein.positioned_mol)
        score_file = os.path.join(self.work_path, self.long_name, self.long_name + '.minimised.json')
        with open(score_file, 'w') as w:
            json.dump({'Energy': energy,
                       'mRMSD': mrsmd.mrmsd,
                       'RMSDs': mrsmd.rmsds}, w)
        Chem.MolToMolFile(ligand, lig_file)
        self._log_warnings()

    # =================== Other ========================================================================================

    @classmethod
    def add_constraint_to_warhead(cls, name:str, constraint:str):
        """
        Add a constraint (multiline is fine) to a warhead definition.
        This will be added and run by Egor's minimiser.

        :param name:
        :param constraint:
        :return: None
        """
        for war_def in cls.warhead_definitions:
            if war_def['name'] == name:
                war_def['constraint'] = constraint
                break
        else:
            raise ValueError(f'{name} not found in warhead_definitions.')


    @classmethod
    def closest_hit(cls, pdb_filenames: List[str],
                 target_resi: int,
                 target_chain: str,
                 target_atomname: str,
                 ligand_resn='LIG') -> str:
        """
        This classmethod helps choose which pdb based on which is closer to a given atom.

        :param pdb_filenames:
        :param target_resi:
        :param target_chain:
        :param target_atomname:
        :param ligand_resn:
        :return:
        """
        best_d = 99999
        best_hit = -1
        with pymol2.PyMOL() as pymol:
            for hit in pdb_filenames:
                pymol.cmd.load(hit)
                d = min([pymol.cmd.distance(f'chain {target_chain} and resi {target_resi} and name {target_atomname}',
                                            f'resn {ligand_resn} and name {atom.name}') for atom in
                         pymol.cmd.get_model(f'resn {ligand_resn}').atom])
                if d < best_d:
                    best_hit = hit
                    best_d = d
                pymol.cmd.delete('*')
        return best_hit

    @classmethod
    def make_covalent(cls, smiles:str, warhead_name:Optional[str]=None) -> Union[str, None]:
        """
        Convert a unreacted warhead to a reacted one in the SMILES

        :param smiles: unreacted SMILES
        :param warhead_name: name in the definitions. If unspecified it will try and guess (less preferrable)
        :return: SMILES
        """
        mol = Chem.MolFromSmiles(smiles)
        if warhead_name:
            war_defs = [wd for wd in cls.warhead_definitions if wd['name'] == warhead_name]
        else:
            war_defs = cls.warhead_definitions
        for war_def in war_defs:
            ncv = Chem.MolFromSmiles(war_def['noncovalent'])
            cv = Chem.MolFromSmiles(war_def['covalent'])
            if mol.HasSubstructMatch(ncv):
                x = Chem.ReplaceSubstructs(mol, ncv, cv, replacementConnectionPoint=0)[0]
                return Chem.MolToSmiles(x)
        else:
            return None

    @classmethod
    def make_all_warhead_combinations(cls, smiles:str, warhead_name:str) -> Union[str, None]:
        """
        Convert a unreacted warhead to a reacted one in the SMILES

        :param smiles: unreacted SMILES
        :param warhead_name: name in the definitions
        :return: SMILES
        """
        mol = Chem.MolFromSmiles(smiles)
        war_def = [wd for wd in cls.warhead_definitions if wd['name'] == warhead_name][0]
        ncv = Chem.MolFromSmiles(war_def['noncovalent'])
        if mol.HasSubstructMatch(ncv):
            combinations = {}
            for wd in cls.warhead_definitions:
                x = Chem.ReplaceSubstructs(mol, ncv, Chem.MolFromSmiles(wd['covalent']), replacementConnectionPoint=0)
                combinations[wd['name']+'_covalent'] = Chem.MolToSmiles(x[0])
                x = Chem.ReplaceSubstructs(mol, ncv, Chem.MolFromSmiles(wd['noncovalent']), replacementConnectionPoint=0)
                combinations[wd['name'] + '_noncovalent'] = Chem.MolToSmiles(x[0])
            return combinations
        else:
            return None

    # =================== Logging ======================================================================================

    @classmethod
    def enable_stdout(cls, level=logging.INFO) -> None:
        """
        The ``cls.journal`` is output to the terminal.

        :param level: logging level
        :return: None
        """
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(level)
        handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
        cls.journal.addHandler(handler)
        # logging.getLogger('py.warnings').addHandler(handler)

    @classmethod
    def enable_logfile(cls, filename='reanimation.log', level=logging.INFO) -> None:
        """
        The journal is output to a file.

        :param filename: file to write.
        :param level: logging level
        :return: None
        """
        handler = logging.FileHandler(filename)
        handler.setLevel(level)
        handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
        cls.journal.addHandler(handler)
        # logging.getLogger('py.warnings').addHandler(handler)

    @classmethod
    def slack_me(cls, msg: str) -> bool:
        """
        Send message to a slack webhook

        :param msg: Can be dirty and unicode-y.
        :return: did it work?
        :rtype: bool
        """
        # sanitise.
        msg = unicodedata.normalize('NFKD', msg).encode('ascii', 'ignore').decode('ascii')
        msg = re.sub('[^\w\s\-.,;?!@#()\[\]]', '', msg)
        r = requests.post(url=os.environ['SLACK_WEBHOOK'],
                          headers={'Content-type': 'application/json'},
                          data=f"{{'text': '{msg}'}}")
        if r.status_code == 200 and r.content == b'ok':
            return True
        else:
            return False

    # =================== Laboratory ===================================================================================

    @classmethod
    def laboratory(cls, entries: List[dict], cores: int = 1):
        pass



