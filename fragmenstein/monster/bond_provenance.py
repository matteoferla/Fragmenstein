from __future__ import annotations
from enum import Enum
from rdkit import Chem
from typing import List

class BondProvenance(Enum):
    """
    Where does the bond come from. This is used to keep names consistent...
    For now original is used in places. The others are interchangeable TBH.
    """
    ORIGINAL = 1 # present in original mols
    MAIN_NOVEL = 2 # joined (zero-atom linking). Closest
    OTHER_NOVEL = 3 # joined (zero-atom linking). Others
    LINKER = 4
    UNASSIGNED = 5 # likely an error

    @classmethod
    def set_all_bonds(cls, mol: Chem.Mol, provenance_name:str) -> None:
        """
        Sets the provenance of all bonds in mol to a category, which is a string from the provenance

        :param mol:
        :param provenance_name: A string original | main_novel " other_novel | linker
        :return:
        """
        for bond in mol.GetBonds():
            cls.set_bond(bond, provenance_name)


    @classmethod
    def get_bond(cls, bond: Chem.Bond) -> BondProvenance:
        if not bond.HasProp('_Provenance'):
            return cls.UNASSIGNED
        else:
            p = bond.GetIntProp('_Provenance')
            return cls(p)

    @classmethod
    def get_bonds(cls, bonds: List[Chem.Bond]) -> List[BondProvenance]:
        return [cls.get_bond(bond) for bond in bonds]

    @classmethod
    def set_bond(cls, bond: Chem.Bond, provenance_name: str) -> None:
        bond.SetIntProp('_Provenance', cls[provenance_name.upper()].value)

    @classmethod
    def set_bonds(cls, bonds: List[Chem.Bond], provenance_name: str) -> None:
        for bond in bonds:
            cls.set_bond(bond, provenance_name)

    @classmethod
    def has_bond(cls, bond: Chem.Bond) -> bool:
        return bond.HasProp('_Provenance')

    @classmethod
    def copy_bond(cls, donor: Chem.Bond, acceptor: Chem.Bond) -> None:
        p = cls.get_bond(donor)
        acceptor.SetIntProp('_Provenance', p.value)



