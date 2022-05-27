from typing import (Optional, Union, Tuple)

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
        plane = plane.lower()
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

    def ring_on_plane(self, mol_idx: int, ring_idx: int = 0, plane: str = 'yx', centroid_to_origin: bool = True):
        """
        Place the ring flat on a plane.
        The maths is shoddy and as a result there may be a bit of offset.

        centroid_to_origin=False may not fully keep the centroid in place.
        """
        # move to centroid for easy maths
        original_centroid = self.get_centroid_of_ring(ring_idx=ring_idx, mol_idx=mol_idx)
        self.atom_to_origin(mol_idx=mol_idx, atom_idx=original_centroid)  # first atom to origin
        for i in self.get_mol(mol_idx).GetRingInfo().AtomRings()[ring_idx]:
            self.atom_on_plane(mol_idx=mol_idx,
                               atom_idx=i,
                               plane=[plane, plane[::-1]][i % 2]  # I am not sure this is needed.
                               )
            centroid = self.get_centroid_of_ring(ring_idx=ring_idx, mol_idx=mol_idx)
            self.atom_to_origin(mol_idx=mol_idx, atom_idx=centroid)  # first atom to origin
        if not centroid_to_origin:
            self.translate_by_point(mol_idx=mol_idx, point=original_centroid, scale=1)

    def flatten_trio(self, mol_idx: int, atom_idcs: Tuple[int, int, int],
                     primary_axis: str = 'x',
                     secondary_axis: str = 'y',
                     first_to_origin: bool = True,
                     ):
        """
        Give three atom indices place the first and third on the primary axis and
        the second on the plane with a second axis.

        Say there's a propane and the axes are 1ary=x and 2ary=y then this method
        will place the triangle along the x axis, flat on the xy plane, with no z elevation.
        """
        original_origin = self.get_point(mol_idx=mol_idx, atom_idx=atom_idcs[0])
        self.atom_to_origin(mol_idx=mol_idx, atom_idx=atom_idcs[0])  # first atom to origin
        self.atom_on_plane(mol_idx=mol_idx, atom_idx=atom_idcs[2],
                           plane=primary_axis + secondary_axis)  # atom 2 flat (z=0)
        self.atom_on_plane(mol_idx=mol_idx, atom_idx=atom_idcs[1],
                           plane=secondary_axis + primary_axis)  # atom 1 flat (z=0)
        self.atom_on_axis(mol_idx=mol_idx, atom_idx=atom_idcs[2], axis=primary_axis)  # atom 2 on x axis
        if not first_to_origin:
            self.translate_by_point(mol_idx=mol_idx, point=original_origin, scale=1)

    # ----- quicker coord getters --------------------------------------------------

    def print_coords(self, mol_idx: int = 0, no_hydrogens: bool = True):
        for atom in self.get_mol(mol_idx).GetAtoms():  #: Chem.Atom  # noqa
            if atom.GetAtomicNum() == 1 and no_hydrogens:
                continue  # no hydrogens
            atom_idx: int = atom.GetIdx()
            coords = self.get_point(atom_idx, mol_idx)
            print(f'atom {atom_idx: >2}: x={coords.x:.1f} y={coords.y:.1f} z={coords.z:.1f}')  # legit print

    def get_centroid_of_atoms(self,
                              *idx_or_points: Union[int, Geometry.Point3D],
                              mol_idx: int = -1) -> Geometry.Point3D:
        """
        centroid of atoms is the euclidean mean of the atoms

        if a whole ring is wanted see ``get_centroid_of_ring``
        """
        points = [self.get_point(ip, mol_idx) for ip in idx_or_points]
        mean_coord: Callable[[int], float] = lambda ai: float(np.mean([p[ai] for p in points]))  # noqa
        return Geometry.Point3D(mean_coord(0), mean_coord(1), mean_coord(2))

    def get_centroid_of_ring(self, mol_idx: int, ring_idx: int = 0) -> Geometry.Point3D:
        """
        centroid of ring is the euclidean mean of the ring's atoms

        if only some atoms are wanted see ``get_centroid_of_atoms``
        """
        return self.get_centroid_of_atoms(*self.get_mol(mol_idx).GetRingInfo().AtomRings()[ring_idx], mol_idx=mol_idx)

    def get_ring_radius(self, mol_idx: int, ring_idx: int = 0) -> float:
        """
        Get the distance to the centroid to the first atom in the ring.
        """
        atoms = self.get_mol(mol_idx).GetRingInfo().AtomRings()[ring_idx]
        centroid = self.get_centroid_of_atoms(*atoms, mol_idx=mol_idx)
        return atoms[0].Distance(centroid)
