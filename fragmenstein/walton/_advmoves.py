from typing import (Optional, Union)

import numpy as np
from rdkit import Geometry

from ._movements import WaltonMove


class WaltonAdvMove(WaltonMove):

    def atom_on_plane(self, mol_idx: int, atom_idx: Union[int, Geometry.Point3D], plane: str):
        """
        Rotate the mol at ``mol_idx`` (``int``) in the attribite ``mols`` (``List[Chem.Mol]``)
        so that atom at ``atom_idx`` is on a plane (eg. 'xy', 'yz' or 'xz').

        It does so by getting the angle subtended by the position on the axis
        that is not on that plane and rotating it away.

        Say, a point on the plane of the axes x and y have a z position of zero.
        A point that is not can be rotated in the opposite direction of
        the angle between the z position (the opposite leg)
        and either the x or y position (the adjecent leg) around the third axis (y or x),
        thus removing the z position.
        The angle in radians is the inverse tangent of
        the ratio of the opposite over the adjencent.
        """
        coords: Geometry.Point3D = self.get_point(atom_idx, mol_idx)
        assert plane in ['xy', 'xz', 'yx', 'yz', 'zx', 'zy']
        axis = plane[1]
        annulled = (set('xyz') - set(plane)).pop()
        opposite = getattr(coords, annulled)
        adjencent = getattr(coords, plane[0])
        if adjencent != 0:  # normal
            theta = -np.arctan(opposite / adjencent)
            self.rotate(mol_idx=mol_idx,
                        theta=theta,
                        axis=axis,
                        degrees=False)
        elif getattr(coords, plane[1]) != 0:  # one is already flat
            self.atom_on_plane(mol_idx, coords, plane[::-1])
        else:  # this is likely the origin
            pass

    def atom_on_axis(self, mol_idx: int, atom_idx: Union[int, Geometry.Point3D], axis: str):
        """

        Rotate the mol at ``mol_idx`` (``int``) in the attribite ``mols`` (``List[Chem.Mol]``)
        so that atom at ``atom_idx`` is on a plane (eg. 'xy', 'yz' or 'xz').
        This calls ``atom_on_plane`` as an axis in on two axis planes
        """
        axis_on_planes = {'x': ('xy', 'xz'),
                          'y': ('xy', 'yz'),
                          'z': ('xz', 'yz'),
                          }
        assert axis in axis_on_planes, f'Axis {axis} not in {axis_on_planes}'
        for plane in axis_on_planes[axis]:
            self.atom_on_plane(mol_idx, atom_idx, plane)

    def atom_to_origin(self, mol_idx: int, atom_idx: Union[int, Geometry.Point3D]):
        """
        Translate to origin.
        """
        coords: Geometry.Point3D = self.get_point(atom_idx, mol_idx)
        self.translate_by_point(mol_idx=mol_idx, point=coords, scale=-1)

    def translate_by_point(self, mol_idx: int, point: Geometry.Point3D, scale: float = +1):
        """
        Translate by the vector point (``rdkit.Geometry.Point3D``).
        Setting scale to -1 will subtract the point as seen
        """
        self.translate(mol_idx=mol_idx,
                       x=scale * point.x,
                       y=scale * point.y,
                       z=scale * point.z)

    def translate_parallel(self,
                           mol_idx: int,
                           distance: float,
                           base_atom_idx: Union[int, Geometry.Point3D],
                           pointer_atom_idx: Union[int, Geometry.Point3D],
                           ref_mol_idx: Optional[int] = None):
        if ref_mol_idx is None:
            ref_mol_idx = mol_idx
        from_coord: Geometry.Point3D = self.get_point(base_atom_idx, ref_mol_idx)
        to_coord: Geometry.Point3D = self.get_point(pointer_atom_idx, ref_mol_idx)
        # This is a unitary vector
        direction: Geometry.Point3D = from_coord.DirectionVector(to_coord)
        self.translate_by_point(mol_idx=mol_idx, point=direction, scale=distance)

