from IPython.display import display
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, Draw, rdMolAlign
from rdkit import Geometry
from typing import (Any, Callable, ClassVar, ForwardRef, Generic, Optional, Tuple, Type, TypeVar, Union, AbstractSet,
                    ByteString, Container, ContextManager, Hashable, ItemsView, Iterable, Iterator, KeysView, Mapping,
                    MappingView, MutableMapping, MutableSequence, MutableSet, Sequence, Sized, ValuesView, Awaitable,
                    AsyncIterator, AsyncIterable, Coroutine, Collection, AsyncGenerator, AsyncContextManager,
                    Reversible, SupportsAbs, SupportsBytes, SupportsComplex, SupportsFloat, SupportsInt, SupportsRound,
                    ChainMap, Counter, Deque, Dict, DefaultDict, List, OrderedDict, Set, FrozenSet, NamedTuple,
                    Generator, AnyStr, cast, get_type_hints, NewType, no_type_check, no_type_check_decorator, NoReturn,
                    overload, Text, TYPE_CHECKING)
import nglview as nv
from io import StringIO
from fragmenstein import Monster
from fragmenstein.branding import divergent_colors

import numpy as np

from ._base import WaltonBase


class WaltonMove(WaltonBase):

    # get coordinates
    def get_point(self, atom_idx: Union[int, Geometry.Point3D], mol_idx: int = -1) -> Geometry.Point3D:
        """
        Given atom index and mol index return a point.
        If a point itself is given, return that.
        """
        if isinstance(atom_idx, Geometry.Point3D):
            return atom_idx
        elif mol_idx == -1:
            raise ValueError('If providing an atom index as a number please provide a mol_idx')
        elif isinstance(atom_idx, int):
            return self.get_mol(mol_idx).GetConformer(0).GetAtomPosition(atom_idx)
        else:
            raise TypeError(atom_idx)

    # ----------- superpose ------------------------------------------------------------------

    def superpose_by_map(self, maps: Dict[Tuple[int, int], Dict[int, int]]):
        """
        The map is a dict of a tuple of two indices: the moving one and the fixed mol.
        The items are dictionary of atom indices from the former and the correspondance in the latter.
        """
        for (moved_idx, fixed_idx), atom_map in maps.items():
            Chem.rdMolAlign.AlignMol(self.mols[moved_idx],
                                     self.mols[fixed_idx],
                                     atomMap=list(atom_map.items()))
        self.superposeed = True

    def superpose_by_mcs(self, fixed_mol_idx=-1, **mcs_settings):
        """
        Aligns by MCS. The fix molecule is the last molecule unless specified.
        """
        fixed = self.get_mol(fixed_mol_idx)
        commons = []
        for moved in self.mols[:1]:
            res: rdFMCS.MCSResult = rdFMCS.FindMCS([moved, fixed], **mcs_settings)
            common = Chem.MolFromSmarts(res.smartsString)
            Chem.SanitizeMol(common)
            mcs_map = list(zip(moved.GetSubstructMatch(common),
                               fixed.GetSubstructMatch(common)
                               )
                           )
            Chem.rdMolAlign.AlignMol(moved, fixed, atomMap=mcs_map)
            commons.append(common)
        self.superposeed = True
        return commons

    # ----- transforms ----------------------------------------

    def transform(self, mol_idx: int, transform_matrix: np.array):
        """
        Applies an affine transform 4x4 matrix to mol_idx.
        Use ``translate`` or ``rotate``.
        """
        AllChem.TransformConformer(self.get_mol(mol_idx).GetConformer(0),
                                   transform_matrix.astype(np.double)
                                   )

    def translate(self, mol_idx: int, x: float, y: float, z: float):
        """
        Move the molecule at index ``mol_idx`` of ``.mols``.
        """
        translation = np.array([[1, 0, 0, x],
                                [0, 1, 0, y],
                                [0, 0, 1, z],
                                [0, 0, 0, 1]], dtype=np.double)
        self.transform(mol_idx, translation)

    def rotate(self, mol_idx: int, theta: float, axis: str, degrees: bool = True):
        """
        Rotate a molecule around a give angle on an axis.
        This is coded terribly.

        To pivot around a given atom  has to translate the pivot point to the origin, rotate and translate back again.
        If this is gibberish, please look up affine transform 4x4 matrices.

        This is not done by Euler angle because its a one off thing.
        """
        correct_unit = lambda d: np.deg2rad(d) if degrees else d  # noqa No lambda PEP8, ha!
        theta_rad = correct_unit(theta)
        # row, column
        ori_rotation = np.zeros((4, 4)) + np.diag(np.ones(4), 0)
        if axis == 'x':  # roll
            first_row = 1
            second_row = 2
        elif axis == 'y':  # pitch
            first_row = 0
            second_row = 2
        elif axis == 'z':  # yaw
            first_row = 0
            second_row = 1
        else:
            raise ValueError(f'Axis {axis} is not x,y or z')
        # for sanity's sake
        first_col = first_row
        second_col = second_row
        ori_rotation[first_row, first_col] = np.cos(theta_rad)
        ori_rotation[first_row, second_col] = -np.sin(theta_rad)
        ori_rotation[second_row, first_col] = np.sin(theta_rad)
        ori_rotation[second_row, second_col] = np.cos(theta_rad)
        self.transform(mol_idx, ori_rotation)

    def experimental_rotate(self, mol_idx: int, angle_x: float, angle_y: float, angle_z: float, degrees: bool = True):
        """
        THIS HAS A BUG.

        Rotate a molecule around a give set of angles,
        in degrees pivoting on the origin.
        These are Tait-Bryan angles.
        There is a bug in here. Use rotate instead.

        To pivot around a given atom  has to translate the pivot point to the origin, rotate and translate back again.
        If this is gibberish, please look up affine transform 4x4 matrices.
        """

        correct_unit = lambda d: np.deg2rad(d) if degrees else d  # noqa
        alpha = correct_unit(angle_x)
        beta = correct_unit(angle_y)
        gamma = correct_unit(angle_z)
        ori_rotation = np.array([[np.cos(beta) * np.cos(gamma),
                                  np.sin(alpha) * np.sin(beta) * np.cos(gamma) - np.cos(alpha) * np.sin(gamma),
                                  np.cos(alpha) * np.sin(beta) * np.cos(gamma) + np.sin(alpha) * np.sin(gamma),
                                  0],
                                 [np.cos(beta) * np.sin(gamma),
                                  np.sin(alpha) * np.sin(beta) * np.sin(gamma) + np.cos(alpha) * np.cos(gamma),
                                  np.cos(alpha) * np.sin(beta) * np.sin(gamma) - np.sin(alpha) * np.cos(gamma),
                                  0],
                                 [-np.sin(beta),
                                  np.sin(alpha) * np.cos(beta),
                                  np.cos(alpha) * np.cos(beta),
                                  0],
                                 [0, 0, 0, 1]], dtype=np.double)
        self.transform(mol_idx, ori_rotation)
