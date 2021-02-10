"""
Keep a copy of the mol.
"""
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
        self.modifications[label] = copy

    def keep_copies(self, mols: List[Chem.Mol], label=None):
        for i, mol in enumerate(mols):
            copy = Chem.Mol(mol)
            if label is None:
                this_label = f'Mol#{len(self.modifications)}'
            else:
                this_label = f'{label}#{i}'
            copy.SetProp('_Name', this_label)
            self.modifications[this_label] = copy
