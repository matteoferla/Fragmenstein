from ..victor import Victor
from .fritz import Fritz
from ..error import FragmensteinError
from rdkit import Chem
from typing import List, Optional, Dict, Union, Callable
import os, json

class OpenVictor(Victor):
    def _calculate_combination_thermo(self):
        if self.is_covalent:
            raise NotImplementedError('OpenVictor does not support covalent ligands')
        self.journal.debug(f'{self.long_name} - Starting system setup')
        self.mol = Chem.Mol(self.monster.positioned_mol)
        self.fritz = Fritz(positioned_mol=self._get_preminimized_undummied_monster(),
                           pdb_block=self.apo_pdbblock,
                           resn=self.ligand_resn,
                           restrained_atomnames=self._get_restraining_atom_indices(),
                           restraint_k=self.settings.get('restraint_k', 1000),
                           mobile_radius=self.settings.get('mobile_radius', 8.0),
                           )
        self.unminimized_pdbblock = self.fritz.to_pdbblock()
        self._data: Dict = self.fritz.reanimate()
        self._data['origins'] = self.monster.origin_from_mol()
        self.minimized_pdbblock = self.fritz.to_pdbblock()
        self.minimized_mol = self.fritz.to_mol()
        self.checkpoint()

    def _calculate_placement_thermo(self):
        if self.is_covalent:
            raise NotImplementedError('OpenVictor does not support covalent ligands')
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

    def checkpoint(self):
        if not self.settings.get('save_outputs', False):
            return
        self.journal.debug(f'{self.long_name} - saving data to disk at {self.work_path}')
        self.make_output_folder()
        with open(os.path.join(self.work_path, self.long_name + '.holo_unminimised.pdb'), 'w') as w:
            w.write(self.unminimized_pdbblock)
        with open(os.path.join(self.work_path, self.long_name + '.holo_minimised.pdb'), 'w') as w:
            w.write(self.minimized_pdbblock)
        for hit in self.hits:
            Chem.MolToMolFile(hit, os.path.join(self.work_path, hit.GetProp('_Name') + '.mol'))
        Chem.MolToMolFile(self.monster.positioned_mol,
                          os.path.join(self.work_path, self.long_name + '.positioned.mol'))
        Chem.MolToMolFile(self.minimized_mol,
                          os.path.join(self.work_path, self.long_name + '.minimised.mol'))
        with open(os.path.join(self.work_path, self.long_name + '.minimised.json'), 'w') as w:
            json.dump(self._data, w)

