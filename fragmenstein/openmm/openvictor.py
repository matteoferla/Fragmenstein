from ..victor import Victor
from .fritz import Fritz
from ..error import FragmensteinError
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List, Optional, Dict, Union, Callable
import os, json, time

class OpenVictor(Victor):
    uses_pyrosetta = False

    def _calculate_combination_thermo(self):
        return self._calculate_thermo()

    def _calculate_placement_thermo(self):
        return self._calculate_thermo()

    def _calculate_thermo(self):
        if self.is_covalent:
            raise NotImplementedError('OpenVictor does not support covalent ligands')
        self.journal.debug(f'{self.long_name} - Starting system setup')
        self.mol = AllChem.AddHs(self.monster.positioned_mol, addCoords=True)
        restraint_k = self.settings['mm_restraint_k'] # default 1000.
        tolerance = self.settings['mm_tolerance']  # 10 * mmu.kilocalorie_per_mole / (mmu.nano * mmu.meter)
        maxIterations = self.settings['mm_max_iterations'] # 0 is infinite
        self.fritz = Fritz(prepped_mol=self.preminimized_undummied_mol,
                           pdb_block=self.apo_pdbblock,
                           resn=self.ligand_resn,
                           restraining_atom_indices=self._get_restraining_atom_indices(),
                           restraint_k=restraint_k,
                           mobile_radius=self.settings['mm_mobile_radius'],
                           )
        self.unminimized_pdbblock = self.fritz.to_pdbblock()
        self._data: Dict = self.fritz.reanimate(tolerance=tolerance, maxIterations=maxIterations)
        while not self.quick_reanimation and self._data['binding_dG'] > 0 * self.fritz.molar_energy_unit:
            restraint_k = 2.
            self.fritz.alter_restraint(restraint_k / 3.)
            self._data: Dict = self.fritz.reanimate(tolerance=tolerance, maxIterations=maxIterations)
            self.journal.debug(f'{self.long_name} - restraints at {restraint_k}')
        self._data['restraint_k'] = restraint_k
        self._data['origins'] = self.monster.origin_from_mol()
        self.energy_score = {k: v.value_in_unit(self.fritz.molar_energy_unit)
                                    if hasattr(v, 'value_in_unit') else v for k, v in self._data.items()
                                                                 if k not in ('minimized_pdb', 'origins')}
        self.energy_score['unit'] = str(self.fritz.molar_energy_unit)
        self.minimized_pdbblock = self.fritz.to_pdbblock()
        self.minimized_mol = self.fritz.to_mol()
        self.mrmsd = self._calculate_rmsd()
        self.checkpoint()
        self.tock = time.time()
        return self.summarize()

    def _process_settings(self):
        pass

    def _get_restraining_atom_indices(self) -> List[int]:
        # Place has '_Stdev': 1.887379141862766e-15, '_Origin': '["x0395.9", "x0434.6"]', '_Max': 0.37818712299601354}
        # Combine has '_ori_i': 100, '_ori_name': 'x0395', '_x': 9.309, '_y': -5.402, '_z': 26.27
        # and occassionally '_ring_i'
        # The problem stemming from the collasing business
        # for now, all atoms are the same.
        restrainables: List[int] = []
        for atom in self.monster.positioned_mol.GetAtoms():
            if atom.GetSymbol() == 'H':
                continue
            elif atom.HasProp('_ori_name') or atom.HasProp('_Origin'):
                restrainables.append( atom.GetIdx() )
            else:
                # no restraint
                pass
        return restrainables

    def checkpoint(self, save_outputs: Optional[bool] = None):
        if save_outputs is True:
            pass
        elif self.settings.get('save_outputs', False):
            return
        self.journal.debug(f'{self.long_name} - saving data to disk at {self.work_path}')
        self.make_output_folder()
        with open(os.path.join(self.work_path, self.long_name, self.long_name + '.holo_unminimised.pdb'), 'w') as w:
            w.write(self.unminimized_pdbblock)
        with open(os.path.join(self.work_path, self.long_name, self.long_name + '.holo_minimised.pdb'), 'w') as w:
            w.write(self.minimized_pdbblock)
        for hit in self.hits:
            Chem.MolToMolFile(hit, os.path.join(self.work_path, self.long_name, hit.GetProp('_Name') + '.mol'))
        Chem.MolToMolFile(self.monster.positioned_mol,
                          os.path.join(self.work_path, self.long_name, self.long_name + '.positioned.mol'))
        Chem.MolToMolFile(self.fritz.prepped_mol,
                          os.path.join(self.work_path, self.long_name, self.long_name + '.prepped.mol'))
        Chem.MolToMolFile(self.minimized_mol,
                          os.path.join(self.work_path, self.long_name, self.long_name + '.minimised.mol'))
        with open(os.path.join(self.work_path, self.long_name, self.long_name + '.minimised.json'), 'w') as w:
            json.dump({**self.energy_score, 'origins': self._data['origins']}, w)

    def summarize(self):
        if self.error_msg:
            if self.monster is None:
                N_constrained_atoms = float('nan')
                N_unconstrained_atoms = float('nan')
            elif self.monster.positioned_mol is None:
                N_constrained_atoms = float('nan')
                N_unconstrained_atoms = float('nan')
            else:
                N_constrained_atoms = self.constrained_atoms
                N_unconstrained_atoms = self.unconstrained_heavy_atoms
            return {'name': self.long_name,
                    'smiles': self.smiles,
                    'error': self.error_msg,
                    'mode': self.merging_mode,
                    '∆∆G': float('nan'),
                    '∆G_bound': float('nan'),
                    '∆G_unbound': float('nan'),
                    'comRMSD': float('nan'),
                    'N_constrained_atoms': N_constrained_atoms,
                    'N_unconstrained_atoms': N_unconstrained_atoms,
                    'runtime': self.tock - self.tick,
                    'regarded': self.monster.matched,
                    'disregarded': self.monster.unmatched
                    }
        else:
            return {'name': self.long_name,
                    'smiles': self.smiles,
                    'error': self.error_msg,
                    'mode': self.merging_mode,
                    '∆∆G': self.energy_score.get('binding_dG', float('nan')),
                    'comRMSD': self.mrmsd.mrmsd,
                    'N_constrained_atoms': self.constrained_atoms,
                    'N_unconstrained_atoms': self.unconstrained_heavy_atoms,
                    'runtime': self.tock - self.tick,
                    'regarded': self.monster.matched,
                    'disregarded': self.monster.unmatched,
                    'origins': self._data['origins'],
                    }

    def _fix_minimized(self, *args, **kwargs):
        raise NotImplementedError('OpenVictor does not support covalent ligands')
