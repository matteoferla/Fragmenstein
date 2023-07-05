from rdkit import Chem
from ..victor import Victor
from ..m_rmsd import mRMSD
import os, json
from rdkit_to_params import Params


class Wictor(Victor):
    """
    This Victor does not call Igor

    ... code-block:: python
        from framgenstein.faux_victors import Wictor
        from fragmenstein.demo import Mac1

        Wictor.capture_rdkit_log()
        # Wictor.enable_stdout()
        Wictor.error_to_catch = ()
        wicky = Wictor(hits=Mac1.get_n_filtered_mols(2),
                       pdb_block=Mac1.get_template(),
                       )
        wicky.combine()
    """

    def _calculate_combination_thermo(self):
        # override igor.
        pass

    def _calculate_placement_thermo(self):
        # override igor.
        pass

    def _checkpoint_charlie(self):
        # making folder.
        self.make_output_folder()
        # no igor!
        self._log_warnings()
        self.journal.debug(f'{self.long_name} - saving pose collage')
        min_file = os.path.join(self.work_path, self.long_name, self.long_name + '.holo_minimised.pdb')
        with open(min_file, 'w') as w:
            w.write(self.minimized_pdbblock)
        self.journal.debug(f'{self.long_name} - saving Gibbs')
        # recover bonds
        lig_file = os.path.join(self.work_path, self.long_name, self.long_name + '.minimised.mol')
        Chem.MolToMolFile(self.minimized_mol, lig_file)
        score_file = os.path.join(self.work_path, self.long_name, self.long_name + '.minimised.json')
        with open(score_file, 'w') as w:
            json.dump({'Energy': self.energy_score,
                       'mRMSD': self.mrmsd.mrmsd,
                       'RMSDs': self.mrmsd.rmsds}, w)
        self._log_warnings()

    def post_monster_step(self):
        # this is a black overridable methods that will be the last thing called
        self.params = Params.from_mol(self.monster.positioned_mol, name=self.ligand_resn, generic=False)
        self.params.NAME = self.ligand_resn  # force it.
        self.params.polish_mol()
        self.params.comments.clear()
        self.params.comments.append('Generated via Fragmenstein')
        mol = Chem.Mol(self.monster.positioned_mol)
        # allow_lax reduces the constraints if it fails.
        successful: bool = self.monster.mmff_minimize(mol, allow_lax=False)
        self.minimized_mol: Chem.Mol = mol
        self.minimized_pdbblock: str = self._plonk_monster_in_structure()
        # The ddG is how strained the molecule is out of the protein... not the drop from binding.
        self.ddG: float = self.monster.MMFF_score(self.minimized_mol, delta=True)
        self.mrmsd: mRMSD = self._calculate_rmsd()
        # save to disc
        self._checkpoint_charlie()

