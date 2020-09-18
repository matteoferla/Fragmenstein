from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D
from typing import Tuple, List, Dict, Optional, Union
import numpy as np
from warnings import warn
import logging

log = logging.getLogger('Fragmenstein')

class _FragmensteinJoinNeighMixin:
    def join_neighboring_mols(self, mol_A: Chem.Mol, mol_B: Chem.Mol):
        """
        Joins two molecules by first calling _find_closest to find closest.
        That method does all the thinking.
        then by calling _join_atoms.

        :param mol_A:
        :param mol_B:
        :return:
        """
        # offset to avoid clashes in bond by history mode
        self.offset(mol_B)
        # get closets atoms
        combo, anchor_A, anchor_B, distance = self._find_closest(mol_A, mol_B)
        mol = self._join_atoms(combo, anchor_A, anchor_B, distance)
        mol.SetProp('_Name', mol_A.GetProp('_Name') + '~' + mol_B.GetProp('_Name'))
        return mol


    def _join_atoms(self, combo: Chem.RWMol, anchor_A: int, anchor_B: int, distance: float):
        """
        extrapolate positions between.

        :param combo:
        :param anchor_A:
        :param anchor_B:
        :param distance:
        :return:
        """
        conf = combo.GetConformer()
        pos_A = conf.GetAtomPosition(anchor_A)
        pos_B = conf.GetAtomPosition(anchor_B)
        n_new = int(round(distance / 1.22) - 1)
        xs = np.linspace(pos_A.x, pos_B.x, n_new + 2)[1:-1]
        ys = np.linspace(pos_A.y, pos_B.y, n_new + 2)[1:-1]
        zs = np.linspace(pos_A.z, pos_B.z, n_new + 2)[1:-1]

        # correct for ring marker atoms
        def is_ring_atom(anchor: int) -> bool:
            atom = combo.GetAtomWithIdx(anchor)
            if atom.HasProp('_ori_i') and atom.GetIntProp('_ori_i') == -1:
                return True
            else:
                return False

        if is_ring_atom(anchor_A):
            distance -= 1.35 + 0.2 # Arbitrary + 0.2 to compensate for the ring not reaching (out of plane).
            n_new -= 1
            xs = xs[1:]
            ys = ys[1:]
            zs = zs[1:]

        if is_ring_atom(anchor_B):
            distance -= 1.35 + 0.2 # Arbitrary + 0.2 to compensate for the ring not reaching  (out of plane).
            n_new -= 1
            xs = xs[:-1]
            ys = ys[:-1]
            zs = zs[:-1]

        # notify that things could be leary.
        if distance < 0:
            log.info(f'Two ring atoms detected to be too close. Joining for now. They will be bonded/fused/spiro afterwards')
        # place new atoms
        if distance < self.joining_cutoff:
            log.debug(f'Molecules will be joined via atoms {anchor_A}+{anchor_B} ({distance} Å) via the addition of {n_new} atoms.')
            previous = anchor_A
            for i in range(n_new):
                idx = combo.AddAtom(Chem.Atom(6))
                new = combo.GetAtomWithIdx(idx)
                new.SetBoolProp('_Novel', True)
                new.SetIntProp('_ori_i', 999)
                conf.SetAtomPosition(idx, Point3D(float(xs[i]), float(ys[i]), float(zs[i])))
                combo.AddBond(idx, previous, Chem.BondType.SINGLE)
                previous = idx
            combo.AddBond(previous, anchor_B, Chem.BondType.SINGLE)
        else:
            msg = f'Atoms {anchor_A}+{anchor_B} are {distance} Å away. Cutoff is {self.joining_cutoff}.'
            log.warning(msg)
            raise ConnectionError(msg)
        return combo.GetMol()

    # === Find closest =================================================================================================
    # dep methods:

    @staticmethod
    def _closest__is_fullbonded(atom):
        if atom.GetIsAromatic() and len(atom.GetNeighbors()) > 2:
            return True
        elif len(atom.GetNeighbors()) > 3:
            return True
        else:
            return False

    @staticmethod
    def _closest__is_ring_atom(atom):
        if atom.GetIsAromatic():
            return True
        elif atom.HasProp('_ori_i') and atom.GetIntProp('_ori_i') == -1:
            if atom.HasProp('_bonds') and 'AROMATIC' in atom.GetProp('_bonds'):
                return True  # technically it could be non-aromatic (ring fusion).
            else:
                return False
        else:
            return False

    @staticmethod
    def _closest__is_warhead_marked(atom):
        return atom.HasProp('_Warhead') and atom.GetBoolProp('_Warhead') == True

    # https://stackoverflow.com/questions/41921255/staticmethod-object-is-not-callable
    closeness_weights =[
                        (_closest__is_warhead_marked.__func__, np.nan),
                        (_closest__is_fullbonded.__func__, 1.0),
                        (_closest__is_ring_atom.__func__, 0.5)
                       ]

    def _find_closest(self, mol_A: Chem.Mol, mol_B: Chem.Mol) -> Tuple[Chem.RWMol, int, int, float]:
        """
        first step in join_neighboring_mols
        This is not to be confused with cls.find_closest

        :param mol_A:
        :param mol_B:
        :return:
        """

        combo = Chem.RWMol(Chem.CombineMols(mol_A, mol_B))
        # ========= distance matrix pre-tweaks.
        distance_matrix = self._get_distance_matrix(combo, mol_A, mol_B)
        penalties = self._get_penalties(combo, distance_matrix.shape)
        # ========= get closest
        pendist_matrix = penalties + distance_matrix
        pendistance = np.nanmin(pendist_matrix)
        if np.isnan(pendistance):
            raise ConnectionError('This is impossible. Previous is absent??')
        else:
            p = np.where(pendist_matrix == pendistance)
            anchor_A = int(p[0][0])
            anchor_B = int(p[1][0])
            distance = distance_matrix[anchor_A, anchor_B]
            penalty = distance_matrix[anchor_A, anchor_B]
            log.debug(f'Connecting {anchor_A} with {anchor_B}, {penalty} penalised distance of {distance}')
            return combo, anchor_A, anchor_B, distance

    def _get_distance_matrix(self, combo: Chem.Mol, A: Union[Chem.Mol, np.ndarray], B: Union[Chem.Mol, np.ndarray]) -> np.ndarray:
        """
        Called by ``_find_closest`` and ``_determine_mergers_novel_ringcore_pair`` in collapse ring (for expansion).

        :param combo:
        :param mol_A:
        :param mol_B:
        :return:
        """
        # TODO move to base once made.
        # input type
        if isinstance(A, Chem.Mol):
            mol_A = A
            A_idxs = np.arange(mol_A.GetNumAtoms())
        else:
            mol_A = None
            A_idxs = np.array(A)
        if isinstance(B, Chem.Mol):
            mol_B = B
            B_idxs = np.arange(mol_B.GetNumAtoms()) + mol_A.GetNumAtoms()
        else:
            mol_B = None
            B_idxs = np.array(B)
        # make matrix
        distance_matrix = Chem.Get3DDistanceMatrix(combo)
        length = combo.GetNumAtoms()
        # nan fill the self values
        self._nan_submatrix(distance_matrix, A_idxs)
        self._nan_submatrix(distance_matrix, B_idxs)
        return distance_matrix

    def _get_penalties(self, combo: Chem.Mol, shape: Tuple[int, int]) -> np.ndarray:
        """
        Called by ``_find_closest``.

        :param combo:
        :param shape:
        :return:
        """
        # penalise
        penalties = np.zeros(shape)
        for fun, weight in self.closeness_weights:
            weigh_bool = np.array([fun(atom) for atom in combo.GetAtoms()])
            penalties[weigh_bool, :] += weight
            penalties[:, weigh_bool] += weight
        return penalties

    def _nan_submatrix(self, matrix, indices):
        """
        Given a square matrix, blank the self-submatrix of the group of indices
        There is probably a better way to do this.

        :param matrix:
        :param indices:
        :return:
        """
        dimension = matrix.shape[0]
        bool_vector = np.zeros((dimension, 1)).astype(bool)
        bool_vector[indices] = True
        bool_matrix = np.tile(bool_vector, (1, dimension))
        logic = np.logical_and(bool_matrix, bool_matrix.transpose())
        matrix[logic] = np.nan

    # ORIGINAL CODE... TO BE DELETED.
        # previous = None
        # while distance == float('inf'):
        #     distance = np.nanmin(tm2)
        #     p = np.where(tm2 == distance)
        #     anchor_A = int(A_idxs[p[0][0]])
        #     anchor_B = int(B_idxs[p[1][0]])
        #     atom_A = combo.GetAtomWithIdx(anchor_A)
        #     atom_B = combo.GetAtomWithIdx(anchor_B)
        #     def outstrike_A(p):
        #         tm2[p[0][0], :] = np.ones(tm2.shape[1]) * float('inf')
        #     def outstrike_B(p):
        #         tm2[:, p[1][0]] = np.ones(tm2.shape[0]) * float('inf')
        #
        #     if distance == float('inf') and previous is not None:
        #         log.debug(f'run out of options for distance bonding, using previous')
        #         d, anchor_A, anchor_B = previous
        #     elif distance == float('inf'):
        #         raise ConnectionError('This is impossible. Previous is absent??')
        #     elif is_warhead_marked(atom_A):
        #         log.debug(f'Atom A ({anchor_A}) is warhead.')
        #         # previous = (d, anchor_A, anchor_B) # forbid.
        #         distance = float('inf')
        #         outstrike_A(p)
        #     elif is_warhead_marked(atom_B):
        #         log.debug(f'Atom B ({anchor_B}) is warhead.')
        #         # previous = (d, anchor_A, anchor_B) # forbid.
        #         distance = float('inf')
        #         outstrike_B(p)
        #     elif is_fullbonded(atom_A):
        #         log.debug(f'Atom A ({anchor_A}) is already at full bond allowance')
        #         previous = (d, anchor_A, anchor_B)
        #         distance = float('inf')
        #         outstrike_A(p)
        #     elif is_fullbonded(atom_B):
        #         log.debug(f'Atom B ({anchor_B}) is already at full bond allowance')
        #         previous = (d, anchor_A, anchor_B)
        #         distance = float('inf')
        #         outstrike_B(p)
        #     elif is_ring_atom(atom_A) and previous is None:
        #         log.debug(f'Atom A ({anchor_A}) is a ring. Don\'t really want to connect to that')
        #         previous = (d, anchor_A, anchor_B)
        #         distance = float('inf')
        #         outstrike_A(p)
        #     elif is_ring_atom(atom_B) and previous is None:
        #         log.debug(f'Atom B ({anchor_B}) is a ring. Don\'t really want to connect to that')
        #         previous = (d, anchor_A, anchor_B)
        #         distance = float('inf')
        #         outstrike_B(p)
        #     elif previous is not None and previous[0] - distance > 0.5:
        #         log.info(f'going with the ring then, the next one is too far.')
        #         d, anchor_A, anchor_B = previous
        #     else:
        #         if self._debug_draw:
        #             print(is_fullbonded(atom_A), is_fullbonded(atom_B), is_ring_atom(atom_A), is_ring_atom(atom_B),
        #                   previous)
        #         pass
