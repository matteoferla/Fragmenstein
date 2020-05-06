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

from ..egor import Egor
from ..core import Fragmenstein
from typing import List, Tuple, Dict, Union, Optional

from rdkit_to_params import Params, Constraints

import os, logging, sys, re, pymol2

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

    def __init__(self,
                 smiles: str,
                 hits: List[Chem.Mol],
                 pdb_filename: str,
                 long_name: str = 'ligand',
                 ligand_resn: str = 'LIG',
                 ligand_resi: Union[int, str] = '1B',
                 covalent_resn: str = 'CYS',  # no other option is accepted.
                 covalent_resi: Optional[Union[int, str]] = None
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
        # making folder.
        self._make_output_folder()
        # make params
        self.journal.debug(f'{self.long_name} - Starting parameterisation')
        self.params = Params.from_smiles(self.smiles, name=ligand_resn, generic=False)
        self.mol = self.params.mol
        # deal with covalent and non covalent separately
        if self.is_covalent:
            self.journal.debug(f'{self.long_name} - is covalent.')
            self.constraint = self._fix_covalent()
            attachment = self._get_attachment_from_pdbblock()
        else:
            self.journal.debug(f'{self.long_name} - is not covalent.')
            self.constraint = None
            attachment = None
        # ***** FRAGMENSTEIN *******
        # make fragmenstein
        self.journal.debug(f'{self.long_name} - Starting fragmenstein')
        self.fragmenstein = Fragmenstein(self.mol, self.hits, attachment=attachment)
        self.unminimised_pdbblock = self._place_fragmenstein()
        # save stuff
        params_file, holo_file, constraint_file = self._checkpoint()
        # ***** EGOR *******
        self.journal.debug(f'{self.long_name} - setting up Egor')
        self.egor = Egor.from_pdbblock(pdbblock=self.unminimised_pdbblock,
                                       params_file=params_file,
                                       constraint_file=constraint_file,
                                       ligand_residue=self.ligand_resi,
                                       key_residues=[self.covalent_resi])
        self.journal.debug(f'{self.long_name} - Egor minimising')
        self.egor.minimise()
        self.journal.debug(f'{self.long_name} - Complete')

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
                return Constraints(smiles=(war_def['covalent'], cov_def['smiles']),
                                   names=[*war_def['covalent_atomnames'], *cov_def['atomnames']],
                                   covalent_res=self.covalent_resi,
                                   target_res=self.ligand_resi)
        else:
            raise ValueError(f'{self.long_name} - Unsure what the warhead is.')

    def _get_attachment_from_pdbblock(self) -> Chem.Mol:
        """
        Yes, I see the madness in using pymol to get an atom for rdkit to make a pose for pyrosetta.
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
        self.journal.debug(f'{self.long_name} - placing fragmenstein')
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(self.apo_pdbblock, 'apo')
            # distort positions
            pymol.cmd.read_pdbstr(Chem.MolToPDBBlock(self.fragmenstein.positioned_mol), 'scaffold')
            pymol.cmd.remove('name R')  # no dummy atoms!
            pymol.cmd.remove('resn UNL')  # no unmatched stuff.
            pymol.cmd.save(f'{self.work_path}/{self.long_name}/{self.long_name}.scaffold.pdb')
            pdbblock = pymol.cmd.get_pdbstr('*')
            pymol.cmd.delete('*')
        return 'LINK         SG  CYS A 145                 CX  LIG B   1     1555   1555  1.8\n' + pdbblock

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

    def _checkpoint(self):
        #  saving params
        self.journal.debug(f'{self.long_name} - saving params')
        params_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params')
        self.params.dump(params_file)
        # saving holo
        self.journal.debug(f'{self.long_name} - saving holo (unmimised)')
        holo_file = os.path.join(self.work_path, self.long_name, self.long_name + '.pdb')
        with open(holo_file, 'w') as w:
            w.write(self.unminimised_pdbblock)
        # saving constraint
        if self.constraint:
            self.journal.debug(f'{self.long_name} - saving constraint')
            constraint_file = os.path.join(self.work_path, self.long_name, self.long_name + '.con')
            self.constraint.dump(constraint_file)
        else:
            constraint_file = ''
        return params_file, holo_file, constraint_file

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
        cls.journal.addHandler(handler)

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
        cls.journal.addHandler(handler)

    # =================== Laboratory ===================================================================================

    @classmethod
    def laboratory(cls, entries: List[dict], cores: int = 1):
        pass
