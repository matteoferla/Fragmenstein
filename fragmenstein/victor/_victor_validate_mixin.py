from ._victor_base_mixin import _VictorBaseMixin
from ..core import Fragmenstein
from ..igor import Igor
from ..m_rmsd import mRSMD
from ..rectifier import Rectifier
from typing import List, Optional, Dict, Union, Callable
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_to_params import Params, Constraints
import time, warnings, os

class _VictorValidateMixin(_VictorBaseMixin):

    @classmethod
    def validate(cls, hits: List[Chem.Mol], smiles: str, pdb_filename: str, ligand_resn: str = 'LIG'):
        """
        Given a pdb file with the protein with ligand and a smiles string, score it against hits.
        Without changing the position of the ligand —bar for protonations.
        The purpose of this is to get a DDG and a RMSD of a given followup for validation purposes.

        :param hits:
        :param smiles:
        :param pdb_filename:
        :param ligand_resn:
        :return:
        """
        # ***** STORE *******
        # entry attributes
        self = cls.__new__(cls)
        self.fragmenstein_merging_mode = 'none_permissive'
        fname = long_name = os.path.split(pdb_filename)[-1]
        long_name = os.path.splitext(fname)[0] + '_validation'
        self.long_name = self.slugify(long_name)
        self.pdb_filename = pdb_filename # non standard!
        self.unminimised_pdbblock = open(pdb_filename).read()
        self.apo_pdbblock = None
        self.hits = hits
        self.ligand_resn = ligand_resn
        pdb = Chem.MolFromPDBBlock(self.unminimised_pdbblock)
        attachment, attachee = self.find_attachment(pdb, ligand_resn)
        if attachment is not None:
            self.is_covalent = True
            if '*' in smiles:
                self.smiles = smiles
            else:
                self.smiles = self.make_covalent(smiles)
        else:
            self.smiles = smiles
            self.is_covalent = False
            attachment, attachee = self.find_closest_to_ligand(pdb, ligand_resn)
        info = attachment.GetPDBResidueInfo()
        self.covalent_resn = info.GetResidueName()
        self.covalent_resi = str(info.GetResidueNumber())+info.GetChainId()
        info = attachee.GetPDBResidueInfo()
        self.ligand_resi = str(info.GetResidueNumber())+info.GetChainId()
        self.pose_fx = None
        # these are calculated
        self.params = None
        self.mol = None
        self.constraint = None
        self.extra_constraint = None
        self.fragmenstein = None
        self.igor = None
        self.minimised_pdbblock = None
        self.minimised_mol = None
        # this is unique to validate
        self.reference_mol = self.extract_mol(name='crystal',
                                              filepath=pdb_filename,
                                              smiles=smiles,
                                              ligand_resn=ligand_resn)
        # buffers etc.
        self._warned = []
        self.energy_score = {'ligand_ref2015': {'total_score': float('nan')},
                             'unbound_ref2015': {'total_score': float('nan')}}
        self.mrmsd = mRSMD.mock()
        self.tick = time.time()
        self.tock = float('inf')
        # analyse
        self._safely_do(execute=self._vanalyse, resolve=self._resolve, reject=self._reject)
        return self

    def _vanalyse(self):
        # THIS IS A COPY PASTE EXCEPT FOR REANIMATE and Params!!
        #self._assert_inputs()
        # ***** PARAMS & CONSTRAINT *******
        self.journal.info(f'{self.long_name} - Starting work')
        self._log_warnings()
        # making folder.
        self._make_output_folder()
        # make params
        self.journal.debug(f'{self.long_name} - Starting parameterisation')
        self.params = Params.from_smiles_w_pdbfile(pdb_file=self.pdb_filename,
                                                   smiles=self.smiles, generic= False,
                                                   name=self.ligand_resn,
                                                   proximityBonding=False)
        self.journal.warning(f'{self.long_name} - CHI HAS BEEN DISABLED')
        self.params.CHI.data = []  # TODO fix chi
        self.mol = self.params.mol
        self._log_warnings()
        # get constraint
        self.constraint = self._get_constraint(self.extra_constraint)
        attachment = self._get_attachment_from_pdbblock() if self.is_covalent else None
        self._log_warnings()
        self.post_params_step()
        # ***** FRAGMENSTEIN *******
        # make fragmenstein
        attachment = self._get_attachment_from_pdbblock() if self.is_covalent else None
        self.journal.debug(f'{self.long_name} - Starting fragmenstein')
        self.fragmenstein = Fragmenstein(mol=self.params.mol, #Chem.MolFromSmiles(self.smiles)
                                         hits=self.hits,
                                         attachment=attachment,
                                         merging_mode=self.fragmenstein_merging_mode,
                                         debug_draw=self.fragmenstein_debug_draw,
                                         average_position=self.fragmenstein_average_position)
        if self.fragmenstein_mmff_minisation:
            self.journal.debug(f'{self.long_name} - pre-minimising fragmenstein (MMFF)')
            self.fragmenstein.mmff_minimise(self.fragmenstein.positioned_mol)
        self.constraint.custom_constraint += self._make_coordinate_constraints()
        self._checkpoint_bravo()
        # save stuff
        params_file, holo_file, constraint_file = self._save_prerequisites()
        self.post_fragmenstein_step()
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
        # DO NOT DO ddG = self.reanimate()
        ddG = self.quick_renamiate() # put igor to work!
        self.minimised_pdbblock = self.igor.pose2str()
        self.post_igor_step()
        self.minimised_mol = self._fix_minimised()
        self.mrmsd = self._calculate_rmsd()
        self.journal.info(f'{self.long_name} - final score: {ddG} kcal/mol {self.mrmsd.mrmsd}.')
        self._checkpoint_charlie()
        #RMSD against self.reference_mol and docked
        # mRSMD.from_other_annotated_mols(self.minimised_mol, self.hits, self.fragmenstein.positioned_mol)
        #dock = self.igor.dock()
        # dock.dump_pdb()
        # self.igor.mol_from_pose(dock)
        # RMSD again
        self.journal.debug(f'{self.long_name} - Completed')

    def quick_renamiate(self):
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

