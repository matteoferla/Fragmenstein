from ._victor_common import _VictorCommon
from ..monster import Monster
from ..igor import Igor
from ..m_rmsd import mRMSD
from typing import List, Optional, Dict, Union, Callable
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_to_params import Params, Constraints
import time, warnings
from ._victor_mode import VictorMinMode
from ..error import DistanceError


# ================== Main entry=====================================================================================



class _VictorCombine(_VictorCommon):

    def combine(self,
                long_name: Optional[str] = None,
                atomnames: Optional[Dict[int, str]] = None,
                warhead_harmonisation: str = 'first',
                joining_cutoff=5.,  # Å
                extra_ligand_constraint: Union[str] = None,
                ):
        """
         Combines the hits without a template.
        If the class attribute ``monster_throw_on_discard`` is True, it will raise an exception if it cannot.

        The cutoff distance is controlled by class attribute ``monster_joining_cutoff``.
        At present this just adds a hydrocarbon chain, no fancy checking for planarity.

        The hits are collapsed, merged, expanded and bonded by proximity.
        In ``(self.monster.expand_ring(..., bonded_as_original=False)`` changing to True, might work, but most likely won't.

        ``warhead_harmonisation`` fixes the warhead in the hits to be homogeneous.

            * ``keep``. Don't do anything
            * ``none``. strip warheads
            * ``first``. Use first warhead
            * warhead name. Use this warhead.

        :param long_name:
        :param atomnames: an optional dictionary that gets used by ``Params.from_smiles``
        :param warhead_harmonisation: keep | strip | first | chloracetimide | nitrile ...
        :param joining_cutoff:
        :param extra_ligand_constraint:
        :return:
        """
        self.joining_cutoff = joining_cutoff
        self.atomnames = atomnames
        self.warhead_harmonisation = warhead_harmonisation
        if long_name is None:
            self.long_name = '-'.join([h.GetProp('_Name') for h in self.hits])
        else:
            self.long_name = self.slugify(long_name)
        self.add_extra_constraint(extra_ligand_constraint)
        # ## Analyse
        self._safely_do(execute=self._calculate_combination, resolve=self._resolve, reject=self._reject)
        return self

    def _harmonize_warhead_combine(self):
        """
        Runs self.harmonize_warheads on the hits, but also determines covalency

        :return:
        """
        self.journal.debug(f'{self.long_name} - harmonising warheads on hits in "{self.warhead_harmonisation}" mode')
        with warnings.catch_warnings(record=True) as self._warned:
            self.hits = self.harmonize_warheads(self.hits, self.warhead_harmonisation, covalent_form=True)
            self._log_warnings()
        # these are calculated
        starhits = any(['*' in Chem.MolToSmiles(h) for h in self.hits])
        if starhits and (self.covalent_resi is None or self.covalent_resn is None):
            raise ValueError(f'{self.long_name} - is covalent but without known covalent residues')
        elif self.warhead_harmonisation == 'strip':
            self.is_covalent = False
        elif starhits:
            self.is_covalent = True
        else:
            self.is_covalent = False


    def _calculate_combination(self):
        """
        called by ``combine`` within ``_safely_do``
        """
        self._calculate_combination_chem()
        self._calculate_combination_thermo()

    def _calculate_combination_chem(self):
        """
        The rdkit part. Monster is used to combine the hits.
        """
        attachment = self._get_attachment_from_pdbblock() if self.is_covalent else None
        self._harmonize_warhead_combine()
        # TODO Does combine not need attachment??
        self.monster.modifications = self.modifications
        self.monster.combine(keep_all=self.monster_throw_on_discard,
                             collapse_rings=True,
                             joining_cutoff=self.joining_cutoff  # Å
                             )
        self.post_monster_step()  # empty overridable
        if self.monster_throw_on_discard and len(self.monster.unmatched):
            raise DistanceError(hits=self.monster.unmatched)
        self.mol = self.monster.positioned_mol
        self.smiles = Chem.MolToSmiles(self.mol)
        # making folder.
        self.make_output_folder()

    def _calculate_combination_thermo(self):
        # paramterise
        self.journal.debug(f'{self.long_name} - Starting parameterisation')
        # ``Params.load_mol`` does a few things: ``from_mol``, the ``polish_mol`` ad then ``convert_mol``.
        # here it is split up.
        self.params = Params.from_mol(self.mol, name=self.ligand_resn, generic=False)
        self.params.NAME = self.ligand_resn  # force it.
        self.params.polish_mol()
        self.params.comments.clear()
        self.params.comments.append('Generated via Fragmenstein')
        # get constraint
        self.constraint = self._get_constraint(self.extra_constraint)
        self.constraint.custom_constraint += self.make_coordinate_constraints_for_combination()
        # _get_constraint will have changed the names in params.mol so the others need changing too!
        # namely  self.params.rename_by_substructure happened.
        self.mol = Chem.Mol(self.params.mol)
        self.monster.positioned_mol = Chem.Mol(self.mol)
        # those lack correct Hs and names and charge!!
        #self.params.add_Hs()  # already done by rectify
        self.params.add_conformer()
        self.params.convert_mol()
        # self.journal.warning(f'{self.long_name} - CHI HAS BEEN DISABLED')
        # self.params.CHI.data = []  # TODO check if chi fix is okay
        self._log_warnings()
        self.post_params_step()  # empty overridable
        self.mmerging_mode = 'full' #TODO: check if this is intended
        self.unminimized_pdbblock = self._plonk_monster_in_structure()
        params_file, holo_file, constraint_file = self._save_prerequisites()
        self.unbound_pose = self.params.test()
        self._checkpoint_alpha()
        self._checkpoint_bravo()
        self.pre_igor_step()  # empty overridable
        self.igor = Igor.from_pdbblock(pdbblock=self.unminimized_pdbblock,
                                       params_file=params_file,
                                       constraint_file=constraint_file,
                                       ligand_residue=self.ligand_resi,
                                       key_residues=[self.covalent_resi])
        # user custom code.
        if self.pose_fx is not None:
            self.journal.debug(f'{self.long_name} - running custom pose mod.')
            self.pose_fx(self.igor.pose)
        else:
            self.pose_mod_step()  # empty overridable
        # storing a roundtrip
        self.unminimized_pdbblock = self.igor.pose2str()
        # minimise until the ddG is negative.
        self.reanimate_n_store()
        self.journal.debug(f'{self.long_name} - Completed')







