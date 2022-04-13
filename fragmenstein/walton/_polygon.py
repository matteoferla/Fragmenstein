import numpy as np
from typing import Optional
from rdkit import Geometry, Chem

class WaltonPolygon:

    @staticmethod
    def create_polygon(n=6, radius=1.5,
                       x0=0., y0=0., z0=0.,
                       star=False,  # technical polygram?
                       bond_type: Chem.BondType = Chem.BondType.SINGLE,
                       starting_mol: Optional[Chem.Mol] = None,
                       name='polygon'):
        """
        Creates a xy flat polygon or star of radius ``radius`` with ``n`` sides.
        """
        if starting_mol:
            mol = Chem.RWMol(starting_mol)
            conf = starting_mol.GetConformer()
        else:
            mol = Chem.RWMol()
            conf = Chem.Conformer()
        center = -1  # just to avoid #noqas
        if star:
            center: int = mol.AddAtom(Chem.Atom(6))
            conf.SetAtomPosition(center, Geometry.Point3D(x0, y0, z0))
        for i in range(n):
            j = mol.AddAtom(Chem.Atom(6))
            if star:
                mol.AddBond(center, j, bond_type)
            elif j == 0:  # no bond yet
                pass
            else:
                mol.AddBond(j - 1, j, bond_type)
            conf.SetAtomPosition(j, Geometry.Point3D(x0 + radius * np.cos(i * 2 * np.pi / n),
                                                     y0 + radius * np.sin(i * 2 * np.pi / n),
                                                     z0)
                                 )
        if not star:
            # bond last atom to first:
            mol.AddBond(j - n + 1, j, bond_type)
        mol.AddConformer(conf)
        mol = mol.GetMol()
        mol.SetProp('_Name', name)
        return mol
