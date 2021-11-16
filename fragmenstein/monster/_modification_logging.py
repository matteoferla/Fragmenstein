########################################################################################################################

__doc__ = \
    """
Keep a copy of the mol.
    """

########################################################################################################################
import json

from rdkit import Chem
from ._base import _MonsterBase
from typing import List

class _MonsterTracker(_MonsterBase):
    """
    _MonsterBase -> _MonsterTracker -> _MonsterCommunal
    """

    def keep_copy(self, mol: Chem.Mol, label=None):
        copy = Chem.Mol(mol)
        if label is None:
            label = f'Mol#{len(self.modifications)}'
        copy.SetProp('_Name', label)
        if label not in self.modifications:
            self.modifications[label] = copy
        else:
            label += '_'
            self.keep_copy(mol, label)


    def keep_copies(self, mols: List[Chem.Mol], label=None):
        for i, mol in enumerate(mols):
            copy = Chem.Mol(mol)
            if label is None:
                this_label = f'Mol#{len(self.modifications)}'
            else:
                this_label = f'{label}#{i}'
            self.keep_copy(mol, this_label)

    def _add_atom_map_asProp(self, mol, atom_map):
        mol.SetProp("atom_map", json.dumps(atom_map))

    def get_atom_map_fromProp(self, mol):
        if mol.HasProp("atom_map"):
            atom_map = json.loads(mol.GetProp("atom_map"))
            atom_map = { int(key):int(val) for key,val in atom_map.items()}
            return atom_map
        else:
            return None
