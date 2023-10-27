from rdkit import Chem
from ..victor import Victor
from ..m_rmsd import mRMSD
import os, json
from rdkit_to_params import Params
from rdkit import Chem
from rdkit.Chem import AllChem


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
    uses_pyrosetta = False

    def _process_settings(self):
        self.journal.debug('Valid settings.py: ff_max_displacement, ff_constraint=1, ff_max_iterations')

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
        if self.settings['ff_use_neighborhood']: # default True
            neighborhood = self.monster.get_neighborhood(self.apo_pdbblock, cutoff=self.settings['ff_neighborhood'])
        else:
            neighborhood = None
        # allow_lax reduces the constraints if it fails
        min_result = self.monster.mmff_minimize(mol,
                                                neighborhood=neighborhood,
                                                ff_max_displacement=float(self.settings['ff_max_displacement']), # def 0
                                                ff_constraint=int(self.settings['ff_constraint']), # def 10
                                                ff_max_iterations=int(self.settings['ff_max_iterations']), # def 200
                                                allow_lax=True)
        self.minimized_mol: Chem.Mol = min_result.mol
        self.minimized_pdbblock: str = self._plonk_monster_in_structure(prepped_mol=self.minimized_mol)
        # The ddG is how strained the molecule is out of the protein... not the drop from binding.
        # recalculating:
        # min_result.ideal is with ff_minimise_ideal True
        ideal: Chem.Mol = self.monster.make_ideal_mol(ff_minimise=bool(self.settings['ff_minimise_ideal']))
        ideal_E: float = ideal.GetDoubleProp('Energy')
        ligand_E: float = self.minimized_mol.GetDoubleProp('Energy')
        # The holo needs recalculating as I don't want the constraints
        AllChem.SanitizeMol(neighborhood)
        hydroneighborhood = Chem.AddHs(neighborhood, addCoords=True)
        holo_E = self.monster.MMFF_score(Chem.CombineMols(self.minimized_mol, hydroneighborhood), delta=False)
        apo_E = self.monster.MMFF_score(hydroneighborhood, delta=False)
        # store data:
        self.energy_score['bound'] = dict(total_score=holo_E, unit='kcal/mol')
        self.energy_score['unbound'] = dict(total_score=apo_E + ideal_E, unit='kcal/mol')
        self.energy_score['apo'] = dict(total_score=apo_E, unit='kcal/mol')
        self.energy_score['ideal'] = dict(total_score=ideal_E, unit='kcal/mol')
        self.energy_score['insitu'] = dict(total_score=ligand_E, unit='kcal/mol')
        self.ddG: float = holo_E - apo_E - ideal_E
        self.mrmsd: mRMSD = self._calculate_rmsd()
        # save to disc
        self._checkpoint_charlie()
