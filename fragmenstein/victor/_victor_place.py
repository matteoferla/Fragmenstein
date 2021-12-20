from typing import *
from ._victor_common import _VictorCommon
from rdkit_to_params import Params
from ..monster import Monster
from ..igor import Igor

class _VictorPlace(_VictorCommon):
    def place(self,
              smiles: str,
              long_name: str = 'ligand',
              merging_mode='none_permissive',
              atomnames: Optional[Dict[int, str]] = None,
              extra_ligand_constraint: Union[str] = None):
        """
        Places a followup (smiles) into the protein based upon the hits.
        Do note that while Monster's place accepts a mol, while place_smiles a smiles
        Victor's place accepts only smiles.

        :param smiles: smiles of followup, optionally covalent (_e.g._ ``*CC(=O)CCC``)
        :param long_name: gets used for filenames so will get corrected
        :param merging_mode:
        :param atomnames: an optional dictionary that gets used by ``Params.from_smiles``
        :param extra_ligand_constraint:
        :return:
        """
        # ## Store
        self.long_name = self.slugify(long_name)
        self.smiles = smiles
        self.atomnames = atomnames
        self.merging_mode = merging_mode
        self.add_extra_constraint(extra_ligand_constraint)
        # ## Analyse
        self._safely_do(execute=self._calculate_placement, resolve=self._resolve, reject=self._reject)
        return self

    def _calculate_placement(self):
        """
        This does all the work

        :return:
        """
        # check they are okay
        self._assert_placement_inputs()
        # ***** PARAMS & CONSTRAINT *******
        self.journal.info(f'{self.long_name} - Starting work')
        self._log_warnings()
        # making folder.
        self.make_output_folder()
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
        self.post_params_step()  # empty overridable
        # ***** FRAGMENSTEIN Monster *******
        # make monster
        self.journal.debug(f'{self.long_name} - Starting fragmenstein')
        # monster_throw_on_discard controls if disconnected.
        self.monster.place(mol=self.mol,
                           attachment=attachment,
                           merging_mode=self.merging_mode)
        self.journal.debug(f'{self.long_name} - Tried {len(self.monster.mol_options)} combinations')
        self.unminimized_pdbblock = self._plonk_monster_in_structure()
        self.constraint.custom_constraint += self.make_coordinate_constraints_for_placement()
        self._checkpoint_bravo()  # saving
        # save stuff
        params_file, holo_file, constraint_file = self._save_prerequisites()
        self.post_monster_step()  # empty overridable
        self.unbound_pose = self.params.test()
        self._checkpoint_alpha()  # saving
        # ***** EGOR *******
        self.journal.debug(f'{self.long_name} - setting up Igor')
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
        if self.quick_reanimation:
            ddG = self.quick_reanimate()
        else:
            ddG = self.reanimate()
        self.ddG = ddG
        self._store_after_reanimation()

    def _assert_placement_inputs(self):
        if '*' in self.smiles and (self.covalent_resi is None or self.covalent_resn is None):
            raise ValueError(f'{self.long_name} - is covalent but without known covalent residues')
            # TODO '*' in self.smiles is bad. user might start with a mol file.
        elif '*' in self.smiles:
            self.is_covalent = True
        else:
            self.is_covalent = False
        assert len(self.ligand_resn) == 3, f'{self.long_name} - {self.ligand_resn} is not 3 char long.'
        assert len(self.hits), f'{self.long_name} - No hits to use to construct'
        assert self.ligand_resn != 'UNL', f'{self.long_name} - It cannot be UNL as it s the unspecified resn in rdkit'
        if self.covalent_resn and len(
                [d for d in self.covalent_definitions if d['residue'] == self.covalent_resn]) == 0:
            raise ValueError(f'{self.long_name} - Unrecognised type {self.covalent_resn}')




