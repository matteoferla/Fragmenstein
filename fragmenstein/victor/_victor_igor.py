from ._victor_store import _VictorStore
from rdkit import Chem
from rdkit.Chem import AllChem
from ..m_rmsd import mRMSD
from typing import Dict, Optional, Callable
from ..extraction_funs import combine_for_bondorder

class _VictorIgor(_VictorStore):

    def _fix_minimized(self, ligand: Optional[Chem.Mol]=None, add_dummy:bool=True) -> Chem.Mol:
        """
        PDBs are terrible for bond order etc. and Rosetta addes these based on atom types
        igor.mol_from_pose cannot use victor.params.mol because it's upstream
        hence this method which adds to it.

        :return:
        """
        if not add_dummy:
            self.journal.debug(f'{self.long_name} - making ligand only')
        else:
            self.journal.debug(f'{self.long_name} - making ligand w/ dummy (if present)')
        if ligand is None:  # normal route
            ligand = self.igor.mol_from_pose(add_dummy=add_dummy)
        new_ligand: Chem.Mol = combine_for_bondorder(self.params.mol, ligand)
        # self.journal.warning(f'{self.long_name} - Rosetta ring closure failed: +infinity kcal/mol penalty')
        # self.energy_score['bound']['total_score'] = float('inf')
        return new_ligand

    def quick_reanimate(self) -> float:
        """
        Correct small deviations from what the forcefield likes.
        Generally flattens buckled rings and that is it.
        Reanimate is normal.

        :return:
        """
        self.igor.coordinate_constraint = 10.
        self.igor.minimize(cycles=5, default_coord_constraint=False)
        self.energy_score = self.calculate_score()
        dG_bound = self.energy_score['bound']['total_score']
        dG_unbound = self.energy_score['unbound']['total_score']
        ddG = dG_bound - dG_unbound
        return ddG

    def reanimate(self) -> float:
        """
        Calls Igor recursively until the ddG is negative or zero.
        igor.minimize does a good job. this is just to get everything as a normal molecule

        :return: ddG (kcal/mol)
        """
        ddG = 999
        self.igor.coordinate_constraint = 0.
        # self.igor.fa_intra_rep = 0.02 # 4x
        # quick unconstrained minimisation to wiggle it out of nasty local minima
        self.igor.minimize(cycles=15, default_coord_constraint=False)
        self.igor.coordinate_constraint = 2
        self.igor.minimize(cycles=5, default_coord_constraint=False)
        self.igor.coordinate_constraint = 1
        while ddG > 0:
            self.journal.debug(f'{self.long_name} - Igor minimising')
            self.igor.minimize(default_coord_constraint=False)
            self.energy_score = self.calculate_score()
            dG_bound = self.energy_score['bound']['total_score']
            dG_unbound = self.energy_score['unbound']['total_score']
            ddG = dG_bound - dG_unbound
            if ddG > 0:
                self.igor.coordinate_constraint /= 2
                self.journal.debug(
                    f'{self.long_name} - coord_constraint lowered: {self.igor.coordinate_constraint}:  {ddG} kcal/mol.')
            if self.igor.coordinate_constraint == 0.:
                self.journal.warning(f'{self.long_name} - failed to minimise without constraints:  {ddG} kcal/mol.')
                break
            elif self.igor.coordinate_constraint < 0.005:
                self.igor.coordinate_constraint = 0.
        self.ddG = ddG
        return ddG

    def reanimate_n_store(self):
        self.reanimate()
        self._store_after_reanimation()

    def _store_after_reanimation(self):
        self.minimized_pdbblock = self.igor.pose2str()
        self.post_igor_step()  # empty overridable
        # extract_mol is a classmethod in _victor_utils: it could be an option but it's inputs differ
        # and there's a few diffent details.
        self.minimized_mol = self._fix_minimized()
        self.mrmsd = self._calculate_rmsd()
        self.journal.info(f'{self.long_name} - final score: {self.ddG} kcal/mol, RMSD: {self.mrmsd.mrmsd}.')
        self._checkpoint_charlie()
        self.journal.debug(f'{self.long_name} - Completed')

    def _calculate_rmsd(self):
        self.journal.debug(f'{self.long_name} - calculating mRMSD')
        return mRMSD.from_other_annotated_mols(self.minimized_mol, self.hits, self.monster.positioned_mol)

    def calculate_score(self):
        return {**self.igor.ligand_score(),
                'unbound': self.igor.detailed_scores(self.unbound_pose, 1)}

    @property
    def constrained_atoms(self) -> int:
        """
        Do note that the whole Origins list contains hydrogens. So do not divided by len!
        :return:
        """
        try:
            conn = sum([o != [] for o in self.monster.origin_from_mol(self.monster.positioned_mol)])
        except KeyboardInterrupt as err:
            raise err
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
        except KeyboardInterrupt as err:
            raise err
        except Exception as err:
            self.journal.warning(f'{self.long_name} - {err.__class__.__name__}: {err}')
            unconn = float('nan')
        return unconn