from __future__ import annotations

########################################################################################################################

__doc__ = \
    """
Victor (after Dr Victor Frankenstein) is a class that uses both Monster (makes blended compounds) and Igor (energy minimises).
This master reanimator keeps a ``.journal`` (logging, class attribute).
And can be called via the class method ``.laboratory`` where he can process multiple compounds at once.

    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.5"
__citation__ = ""

########################################################################################################################

import json
import os
import re
import warnings
import pyrosetta
import time
from typing import List, Union, Optional, Callable, Dict

from rdkit import Chem
from rdkit.Chem import AllChem
# noinspection PyUnresolvedReferences
from rdkit_to_params import Params, Constraints

from ._victor_utils_mixin import _VictorUtilsMixin  # <--- _VictorBaseMixin
from ._victor_validate_mixin import _VictorValidateMixin
from ._victor_automerge_mixin import _VictorAutomergeMixin
from ..monster import Monster
from ..igor import Igor
from ..m_rmsd import mRSMD


class Victor(_VictorUtilsMixin, _VictorValidateMixin, _VictorAutomergeMixin):
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
    ``monster_debug_draw:bool`` and ``monster_merging_mode:str`` are class attributes that control Monster.

    """

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
                 atomnames: Optional[Dict[int, str]] = None
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
        :param atomnames: an optional dictionary that gets used by ``Params.from_smiles``
        """
        # ***** STORE *******
        # entry attributes
        self.long_name = self.slugify(long_name)
        self.smiles = smiles
        with open(pdb_filename) as fh:
            self.apo_pdbblock = fh.read()
        self.hits = hits
        self.ligand_resn = ligand_resn.upper()
        self.ligand_resi = ligand_resi
        self.covalent_resn = covalent_resn.upper()
        self.covalent_resi = covalent_resi
        self.atomnames = atomnames
        self.extra_constraint = extra_constraint
        self.pose_fx = pose_fx
        # these are calculated
        self.is_covalent = None
        self.params = None
        self.mol = None
        self.constraint = None
        self.monster = None
        self.modifications = []  # used by automerger only
        self.unminimised_pdbblock = None
        self.igor = None
        self.minimised_pdbblock = None
        self.minimised_mol = None
        self.reference_mol = None  # filled only for validate
        # buffers etc.
        self._warned = []
        self.energy_score = {'ligand_ref2015': {'total_score': float('nan')},
                             'unbound_ref2015': {'total_score': float('nan')}}
        self.mrmsd = mRSMD.mock()
        self.tick = time.time()
        self.tock = float('inf')
        self.error = ''
        # analyse
        self._safely_do(execute=self._analyse, resolve=self._resolve, reject=self._reject)

    # =================== Init monster methods ============================================================================

    def _safely_do(self,
                   execute: Optional[Callable] = None,
                   resolve: Optional[Callable] = None,
                   reject: Optional[Callable] = None):
        """
        A safety net around the analysis.
        Ought to be a decorator and ought to not use the same names as a JS Promise.
        The class attribute ``error_to_catch`` is by default Exception

        :param execute: what to run (main)
        :param resolve: what to run at the end (regardless of failure)
        :param reject: what to run if ``exceute`` fails
        :return:
        """
        # warnings
        with warnings.catch_warnings(record=True) as self._warned:
            try:
                if execute is not None:
                    execute()
            except self.error_to_catch as err:
                self.error = f'{err.__class__.__name__}: {err}'
                if reject is not None:
                    reject(err)
            finally:
                if resolve is not None:
                    resolve()

    def _resolve(self) -> None:
        """
        This gets called at the end of ``_safely_do``, regardless of the whether it worked or not.
        So the name is a bit misleading.
        :return:
        """
        self.tock = time.time()
        self.journal.info(f'{self.long_name} - Time taken: {self.tock - self.tick}')

    def _reject(self, err) -> None:
        """
        This gets called by ``_safely_do`` on error.

        :param err: the error raised.
        :return:
        """
        self.journal.error(f'{self.long_name} — {err.__class__.__name__}: {err}')

    def _analyse(self) -> None:
        """
        This is the actual monster of the class.

        :return:
        """
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
        self.params = Params.from_smiles(self.smiles, name=self.ligand_resn, generic=False, atomnames=self.atomnames)
        # self.journal.warning(f'{self.long_name} - CHI HAS BEEN DISABLED')
        # self.params.CHI.data = []  # Chi is fixed, but older version. should probably check version
        self.mol = self.params.mol
        self._log_warnings()
        # get constraint
        self.constraint = self._get_constraint(self.extra_constraint)
        attachment = self._get_attachment_from_pdbblock() if self.is_covalent else None
        self._log_warnings()
        self.post_params_step()
        # ***** FRAGMENSTEIN *******
        # make monster
        self.journal.debug(f'{self.long_name} - Starting fragmenstein')
        # monster_throw_on_discard controls if disconnected.
        Monster.throw_on_discard = self.monster_throw_on_discard
        self.monster = Monster(hits=self.hits,
                               debug_draw=self.monster_debug_draw,
                               average_position=self.monster_average_position)
        self.monster.place(mol=self.mol,
                           attachment=attachment,
                           merging_mode=self.monster_merging_mode)
        self.journal.debug(f'{self.long_name} - Tried {len(self.monster.scaffold_options)} combinations')
        self.unminimised_pdbblock = self._place_monster_in_structure()
        self.constraint.custom_constraint += self._make_coordinate_constraints()
        self._checkpoint_bravo()
        # save stuff
        params_file, holo_file, constraint_file = self._save_prerequisites()
        self.post_monster_step()
        self.unbound_pose = self.params.test()
        self._checkpoint_alpha()
        # ***** EGOR *******
        self.journal.debug(f'{self.long_name} - setting up Igor')
        self.igor = Igor.from_pdbblock(pdbblock=self.unminimised_pdbblock,
                                       params_file=params_file,
                                       constraint_file=constraint_file,
                                       ligand_residue=self.ligand_resi,
                                       key_residues=[self.covalent_resi])
        # user custom code.
        if self.pose_fx is not None:
            self.journal.debug(f'{self.long_name} - running custom pose mod.')
            self.pose_fx(self.igor.pose)
        else:
            self.pose_mod_step()
        # storing a roundtrip
        self.unminimised_pdbblock = self.igor.pose2str()
        # minimise until the ddG is negative.
        if self.quick_renanimation:
            ddG = self.quick_reanimate()
        else:
            ddG = self.reanimate()
        self.minimised_pdbblock = self.igor.pose2str()
        self.post_igor_step()
        self.minimised_mol = self._fix_minimised()
        self.mrmsd = self._calculate_rmsd()
        self.journal.info(f'{self.long_name} - final score: {ddG} kcal/mol {self.mrmsd.mrmsd}.')
        self._checkpoint_charlie()
        self.journal.debug(f'{self.long_name} - Completed')

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

    def _make_coordinate_constraints(self):
        """
        See also ``_make_coordinate_constraints_for_unnovels`` in automerge.
        This is the normal function and uses the origin data,
        while the other constrains based on lack of novel attribute.
        """
        lines = []
        origins = self.monster.origin_from_mol(self.monster.positioned_mol)
        std = self.monster.stdev_from_mol(self.monster.positioned_mol)
        mx = self.monster.max_from_mol(self.monster.positioned_mol)
        conf = self.monster.positioned_mol.GetConformer()
        for i in range(self.monster.positioned_mol.GetNumAtoms()):
            if len(origins[i]) > 0:
                atom = self.monster.positioned_mol.GetAtomWithIdx(i)
                if atom.GetSymbol() == '*':
                    continue
                elif atom.GetPDBResidueInfo() is None:
                    self.journal.critical(f'Atom {i} ({atom.GetSymbol()}) has no name!')
                    continue
                pos = conf.GetAtomPosition(i)
                if self.constraint_function_type.upper() == 'HARMONIC':
                    fxn = f'HARMONIC 0 {std[i] + 1}'
                elif self.constraint_function_type.upper() == 'FLAT_HARMONIC':
                    if len(origins[i]) > 1:
                        fxn = f'FLAT_HARMONIC 0 1.0 {mx[i]}'
                    else:
                        fxn = f'HARMONIC 0 1.0'
                elif self.constraint_function_type.upper() == 'BOUNDED':
                    fxn = f'BOUNDED 0 {mx[i]} 1 0.5 TAG'
                else:
                    raise ValueError(f'{self.constraint_function_type} is not HARMONIC or FADE or BOUNDED')
                atomname = atom.GetPDBResidueInfo().GetName()
                lines.append(f'CoordinateConstraint {atomname} {self.ligand_resi} ' + \
                             f'CA {self.covalent_resi} ' + \
                             f'{pos.x} {pos.y} {pos.z} {fxn}\n')
        return ''.join(lines)

    def _get_LINK_record(self):
        if self.is_covalent:
            # get correct chain names.
            l_resi, l_chain = re.match('(\d+)(\D?)', str(self.ligand_resi)).groups()
            if self.covalent_resi:
                p_resi, p_chain = re.match('(\d+)(\D?)', str(self.covalent_resi)).groups()
            else:
                p_resi, p_chain = None, None
            if not p_chain:
                p_chain = 'A'
            if not l_chain:
                l_chain = 'B'
            # get the cx atom name
            cx = self.params.pad_name(self.params.CONNECT[0].atom_name)
            # TODO the SG connection is hardcoded.
            return f'LINK         SG  {self.covalent_resn} {p_chain} {p_resi: >3}                ' + \
                    f'{cx} {self.ligand_resn} {l_chain} {l_resi: >3}     1555   1555  1.8\n'

    def _correct_ligand_info(self, mol: Optional[Chem.Mol]=None) -> Chem.Mol:
        """
        Corrects in place the given mol based on self.ligand_resi.
        If none provided it assumed self.monster.positioned_mol
        """
        if mol is None:
            mol = self.monster.positioned_mol
        l_resi, l_chain = re.match('(\d+)(\D?)', str(self.ligand_resi)).groups() #TODO improve ligand_resi
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            info.SetResidueNumber(int(l_resi))
            info.SetChainId(l_chain)
            info.SetIsHeteroAtom(True)
            info.SetOccupancy(1.)
            info.SetResidueName(self.ligand_resn)
        return mol


    def _place_monster_in_structure(self): # TODO make a new method that does not use Pymol
        """
        Places the molecule in the structure using pymol.
        :return:
        """
        import pymol2
        self._correct_ligand_info()
        mol = AllChem.DeleteSubstructs(self.monster.positioned_mol, Chem.MolFromSmiles('*'))
        if self.monster_mmff_minisation:
            self.journal.debug(f'{self.long_name} - pre-minimising monster (MMFF)')
            self.monster.mmff_minimise(mol)
        self.journal.debug(f'{self.long_name} - placing monster in structure')
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(self.apo_pdbblock, 'apo')
            pos_mol = Chem.MolToPDBBlock(mol)
            pymol.cmd.read_pdbstr(pos_mol, 'scaffold')
            pymol.cmd.remove('name R')  # no dummy atoms!
            for c in self._connected_names:
                pymol.cmd.remove(f'name {c}')  # no conns
            pymol.cmd.remove('resn UNL')  # no unmatched stuff.
            pdbblock = pymol.cmd.get_pdbstr('*')
            pymol.cmd.delete('*')
        return self._get_LINK_record() + pdbblock

    def _fix_minimised(self) -> Chem.Mol:
        """
        PDBs are terrible for bond order etc. and Rosetta addes these based on atom types
        :return:
        """
        self.journal.debug(f'{self.long_name} - making ligand only')
        ligand = self.igor.mol_from_pose()
        template = AllChem.DeleteSubstructs(self.params.mol, Chem.MolFromSmiles('*'))
        return AllChem.AssignBondOrdersFromTemplate(template, ligand)

    def quick_reanimate(self) -> float:
        """
        Correct small deviations from what the forcefield likes. Generally flattens buckled rings and that is it.
        Reanimate is normal.

        :return:
        """
        self.igor.coordinate_constraint = 10.
        self.igor.minimise(cycles=5, default_coord_constraint=False)
        self.energy_score = self.calculate_score()
        dG_bound = self.energy_score['ligand_ref2015']['total_score']
        dG_unbound = self.energy_score['unbound_ref2015']['total_score']
        ddG = dG_bound - dG_unbound
        return ddG

    def reanimate(self) -> float:
        """
        Calls Igor recursively until the ddG is negative or zero.
        igor.minimise does a good job. this is just to get everything as a normal molecule

        :return: ddG (kcal/mol)
        """
        ddG = 999
        self.igor.coordinate_constraint = 0.
        # self.igor.fa_intra_rep = 0.02 # 4x
        # quick unconstrained minimisation to wiggle it out of nasty local minima
        self.igor.minimise(cycles=15, default_coord_constraint=False)
        self.igor.coordinate_constraint = 2
        self.igor.minimise(cycles=5, default_coord_constraint=False)
        self.igor.coordinate_constraint = 1
        while ddG > 0:
            self.journal.debug(f'{self.long_name} - Igor minimising')
            self.igor.minimise(default_coord_constraint=False)
            self.energy_score = self.calculate_score()
            dG_bound = self.energy_score['ligand_ref2015']['total_score']
            dG_unbound = self.energy_score['unbound_ref2015']['total_score']
            ddG = dG_bound - dG_unbound
            if ddG > 0:
                self.igor.coordinate_constraint /= 2
                self.journal.debug(
                    f'{self.long_name} - coord_constraint lowered: {self.igor.coordinate_constraint}:  {ddG} kcal/mol.')
            if self.igor.coordinate_constraint == 0.:
                self.journal.warn(f'{self.long_name} - failed to minimise without constraints:  {ddG} kcal/mol.')
                break
            elif self.igor.coordinate_constraint < 0.005:
                self.igor.coordinate_constraint = 0.
        return ddG

    # =================== Constraint & attachment ======================================================================

    def _get_constraint(self, extra_constraint: Optional[str] = None) -> Union[Constraints, None]:
        # deal with covalent and non covalent separately
        if self.is_covalent:
            self.journal.debug(f'{self.long_name} - is covalent.')
            constraint = self._fix_covalent()
            if extra_constraint:
                constraint.custom_constraint += self.extra_constraint
            return constraint
        else:
            self.journal.debug(f'{self.long_name} - is not covalent.')
            constraint = self._fix_uncovalent()
            if extra_constraint:
                constraint.custom_constraint += self.extra_constraint
            return constraint

    def _fix_uncovalent(self):
        return Constraints.mock()

    def _get_war_def(self):
        for war_def in self.warhead_definitions:
            warhead = Chem.MolFromSmiles(war_def['covalent'])
            if self.mol.HasSubstructMatch(warhead):
                return war_def
        else:
            if self.mol.HasSubstructMatch(Chem.MolFromSmiles('*C')):
                self.journal.warning('Unknown type of covalent')
                return {'name': 'unknown',
                        'covalent': 'C~C*',
                        'covalent_atomnames': ['CY', 'CX', 'CONN1'],
                        'noncovalent': 'C~[C+]',  # clearly not
                        'noncovalent_atomnames': ['CY', 'CX']
                        }
            else:
                raise ValueError(f'{self.long_name} - Unsure what the warhead is {self.smiles}.')

    def _fix_covalent(self):
        self.journal.debug(f'{self.long_name} - fixing for covalent')
        # to make life easier for analysis, CX is the attachment atom, CY is the one before it.
        war_def = self._get_war_def()
        warhead = Chem.MolFromSmiles(war_def['covalent'])
        self.params.rename_by_substructure(warhead, war_def['covalent_atomnames'])
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

    def _get_attachment_from_pdbblock(self) -> Union[None, Chem.Mol]:
        """
        Yes, yes, I see the madness in using pymol to get an atom for rdkit to make a pose for pyrosetta.
        Hence why `find_attachment` will replace it.
        todo `_get_attachment_from_pdbblock` --> `find_attachment`
        """
        import pymol2

        self.journal.debug(f'{self.long_name} - getting attachemnt atom')
        if not self.covalent_resn:
            return None
        else:
            if isinstance(self.covalent_resi, str):
                resi, chain = re.match('(\d+)(\w)', self.covalent_resi).groups()
                resi = int(resi)
            else:
                resi = self.covalent_resi
                chain = None
            with pymol2.PyMOL() as pymol:
                pymol.cmd.read_pdbstr(self.apo_pdbblock, 'prot')
                if self.covalent_resn == 'CYS':
                    name = 'SG'
                else:
                    raise NotImplementedError('only done for cys atm')
                try:
                    if chain is not None:
                        pdb = pymol.cmd.get_pdbstr(f'resi {resi} and name {name} and chain {chain}')
                    else:
                        pdb = pymol.cmd.get_pdbstr(f'resi {resi} and name {name}')
                except:
                    pdb = pymol.cmd.get_pdbstr(f'resi {resi} and name {name}')
                return Chem.MolFromPDBBlock(pdb)

    def _calculate_rmsd(self):
        self.journal.debug(f'{self.long_name} - calculating mRMSD')
        return mRSMD.from_other_annotated_mols(self.minimised_mol, self.hits, self.monster.positioned_mol)

    def calculate_score(self):
        return {**self.igor.ligand_score(),
                'unbound_ref2015': self.igor.detailed_scores(self.unbound_pose, 1)}

    # =================== Other ========================================================================================

    def _log_warnings(self):
        if len(self._warned):
            for w in self._warned:
                self.journal.warning(f'{self.long_name} - {w.message} ({w.category})')
            self._warned.clear()

    # =================== Overridables =================================================================================

    def post_params_step(self):
        """
        This method is intended for make inherited mods easier.
        :return:
        """
        pass

    def post_monster_step(self):
        """
        This method is intended for make inherited mods easier.
        :return:
        """
        pass

    def pose_mod_step(self):
        """
        This method is intended for make inherited mods easier.
        :return:
        """
        pass

    def post_igor_step(self):
        """
        This method is intended for make inherited mods easier.
        :return:
        """
        pass

    # =================== Logging ======================================================================================

    def _save_prerequisites(self):
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
        if self.constraint is not None:
            self.journal.debug(f'{self.long_name} - saving constraint')
            constraint_file = os.path.join(self.work_path, self.long_name, self.long_name + '.con')
            self.constraint.dump(constraint_file)
        else:  # basically impossible.
            constraint_file = ''
        return params_file, holo_file, constraint_file

    def _checkpoint_alpha(self):
        self._log_warnings()
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
        ptest_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params_test.pdb')
        self.unbound_pose.dump_pdb(ptest_file)
        pscore_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params_test.score')
        scorefxn = pyrosetta.get_fa_scorefxn()
        with open(pscore_file, 'w') as w:
            w.write(str(scorefxn(self.unbound_pose)))
        self._log_warnings()

    def _checkpoint_bravo(self):
        self._log_warnings()
        self.journal.debug(f'{self.long_name} - saving mols from monster')
        if self.monster.scaffold is not None:
            scaffold_file = os.path.join(self.work_path, self.long_name, self.long_name + '.scaffold.mol')
            Chem.MolToMolFile(self.monster.scaffold, scaffold_file, kekulize=False)
            if self.monster.scaffold.HasProp('parts'):
                disregard = json.loads(self.monster.scaffold.GetProp('parts'))
                self.journal.info(f'{self.long_name} - disregarded {disregard}')
            else:
                disregard = []
        if self.monster.chimera is not None:
            chimera_file = os.path.join(self.work_path, self.long_name, self.long_name + '.chimera.mol')
            Chem.MolToMolFile(self.monster.chimera, chimera_file, kekulize=False)
        if self.monster.positioned_mol is not None:
            pos_file = os.path.join(self.work_path, self.long_name, self.long_name + '.positioned.mol')
            Chem.MolToMolFile(self.monster.positioned_mol, pos_file, kekulize=False)
        if self.monster.scaffold_options:
            opt_file = os.path.join(self.work_path, self.long_name, self.long_name + '.scaffold_options.sdf')
            writer = Chem.SDWriter(opt_file)
            writer.SetKekulize(False)
            for t in self.monster.scaffold_options:
                writer.write(t)
            writer.close()

        frag_file = os.path.join(self.work_path, self.long_name, self.long_name + '.monster.json')
        data = {'smiles': self.smiles,
                'origin': self.monster.origin_from_mol(self.monster.positioned_mol),
                'stdev': self.monster.stdev_from_mol(self.monster.positioned_mol)}
        if disregard:
            data['disregard'] = disregard
        with open(frag_file, 'w') as w:
            json.dump(data, w)
        self._log_warnings()
        # unminimised_pdbblock will be saved by igor (round trip via pose)

    def _checkpoint_charlie(self):
        self._log_warnings()
        self.journal.debug(f'{self.long_name} - saving pose from igor')
        min_file = os.path.join(self.work_path, self.long_name, self.long_name + '.holo_minimised.pdb')
        self.igor.pose.dump_pdb(min_file)
        self.journal.debug(f'{self.long_name} - calculating Gibbs')
        # recover bonds
        lig_file = os.path.join(self.work_path, self.long_name, self.long_name + '.minimised.mol')
        Chem.MolToMolFile(self.minimised_mol, lig_file)
        score_file = os.path.join(self.work_path, self.long_name, self.long_name + '.minimised.json')
        with open(score_file, 'w') as w:
            json.dump({'Energy': self.energy_score,
                       'mRMSD': self.mrmsd.mrmsd,
                       'RMSDs': self.mrmsd.rmsds}, w)
        self._log_warnings()

    @property
    def constrained_atoms(self) -> int:
        """
        Do note that the whole Origins list contains hydrogens. So do not divided by len!
        :return:
        """
        try:
            conn = sum([o != [] for o in self.monster.origin_from_mol(self.monster.positioned_mol)])
        except Exception as err:
            self.journal.warning(f'{self.long_name} - {err.__class__.__name__}: {err}')
            conn = float('nan')
        return conn

    @property
    def unconstrained_heavy_atoms(self) -> int:
        try:
            origins = self.monster.origin_from_mol(self.monster.positioned_mol)
            unconn = sum([o == [] and atom.GetSymbol() != 'H' for o, atom in
                          zip(origins, self.monster.positioned_mol.GetAtoms())])
        except Exception as err:
            self.journal.warning(f'{self.long_name} - {err.__class__.__name__}: {err}')
            unconn = float('nan')
        return unconn


######### Make params use the same log.

Params.log = Victor.journal
