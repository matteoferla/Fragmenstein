from ._victor_store import _VictorStore
from rdkit import Chem
from rdkit.Chem import AllChem
from ..m_rmsd import mRSMD

class _VictorIgor(_VictorStore):

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
        self.ddG = ddG
        return ddG

    def reanimate_n_store(self):
        self.reanimate()
        self._store_after_reanimation()

    def _store_after_reanimation(self):
        self.minimised_pdbblock = self.igor.pose2str()
        self.post_igor_step()  # empty overridable
        self.minimised_mol = self._fix_minimised()
        self.mrmsd = self._calculate_rmsd()
        self.journal.info(f'{self.long_name} - final score: {self.ddG} kcal/mol, RMSD: {self.mrmsd.mrmsd}.')
        self._checkpoint_charlie()
        self.journal.debug(f'{self.long_name} - Completed')

    def _calculate_rmsd(self):
        self.journal.debug(f'{self.long_name} - calculating mRMSD')
        return mRSMD.from_other_annotated_mols(self.minimised_mol, self.hits, self.monster.positioned_mol)

    def calculate_score(self):
        return {**self.igor.ligand_score(),
                'unbound_ref2015': self.igor.detailed_scores(self.unbound_pose, 1)}

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