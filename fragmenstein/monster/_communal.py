########################################################################################################################
__doc__ = \
    """
This is inherited by all three parents of the place/combine/other group
    """

########################################################################################################################

from typing import Optional, List, Tuple, Union

import numpy as np
from rdkit import Chem

from ._modification_logging import _MonsterTracker
from .bond_provenance import BondProvenance
from ..error import ShoddyCodeError


# _MonsterBase -> _MonsterTracker -> _MonsterCommunal

class _MonsterCommunal(_MonsterTracker):

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
        return atom.HasProp('_Warhead') == 1 and atom.GetBoolProp('_Warhead') is True

    # func: https://stackoverflow.com/questions/41921255/staticmethod-object-is-not-callable
    closeness_weights = [
        (_closest__is_warhead_marked.__func__, np.nan),
        (_closest__is_fullbonded.__func__, 1.0),
        (_closest__is_ring_atom.__func__, 0.5)
        # is_triangle, 2.0
    ]

    def _find_closest(self, mol_A: Chem.Mol, mol_B: Chem.Mol) -> Tuple[Chem.RWMol, int, int, float]:
        """
        first step in join_neighboring_mols
        This is not to be confused with cls.find_closest_to_ligand

        :param mol_A:
        :param mol_B:
        :return:
        """
        combo: Chem.RWMol
        candidates: Tuple[int, int, float]
        combo, candidates = self._find_all_closest(mol_A, mol_B)
        return (combo, *candidates[0])

    def _find_all_closest(self, mol_A: Chem.Mol, mol_B: Chem.Mol) -> Tuple[Chem.RWMol, List[Tuple[int, int, float]]]:
        """
        See _find_closest

        :param mol_A:
        :param mol_B:
        :return:
        """
        combo = Chem.RWMol(Chem.CombineMols(mol_A, mol_B))
        # ========= distance matrix pre-tweaks.
        distance_matrix = self._get_distance_matrix(combo, mol_A, mol_B)
        penalties = self._get_joining_penalties(combo, distance_matrix.shape)
        # ========= get closest
        pendist_matrix = penalties + distance_matrix
        pendistance = float(np.nanmin(pendist_matrix))
        if np.isnan(pendistance):
            raise ShoddyCodeError('This is impossible. Previous is absent??')
        candidates: List[Tuple[int, int, float]] = []

        def get_closest(pendistance: float):
            p = np.where(pendist_matrix == pendistance)
            anchor_A = int(p[0][0])
            anchor_B = int(p[1][0])
            distance = distance_matrix[anchor_A, anchor_B]
            penalty = penalties[anchor_A, anchor_B]
            self.journal.debug(f'Connecting {anchor_A} with {anchor_B} would have a ' +
                               f'{penalty} penalised distance of {distance}')
            return anchor_A, anchor_B, distance

        anchor_A, anchor_B, distance = get_closest(pendistance)
        candidates.append((int(anchor_A), int(anchor_B), distance))
        with np.errstate(invalid='ignore'):
            pendist_matrix[pendist_matrix > 1.] = np.nan
        while pendistance < 1.:
            pendist_matrix[[anchor_A, anchor_B], :] = np.nan
            pendist_matrix[:, [anchor_A, anchor_B]] = np.nan
            pendistance = np.nanmin(pendist_matrix)
            if np.isnan(pendistance):
                break
            else:
                anchor_A, anchor_B, distance = get_closest(pendistance)
                candidates.append((anchor_A, anchor_B, distance))
        return combo, candidates

    def _get_distance_matrix(self,
                             combo: Chem.Mol,
                             A: Union[Chem.Mol, np.ndarray],
                             B: Union[Chem.Mol, np.ndarray]) -> np.ndarray:
        """
        Called by ``_find_closest`` and ``_determine_mergers_novel_ringcore_pair`` in collapse ring (for expansion).
        This is a distance matrix blanked so it is only distances to other fragment

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
        self._nan_fill_submatrix(distance_matrix, list(A_idxs))
        self._nan_fill_submatrix(distance_matrix, list(B_idxs))
        return distance_matrix

    def _nan_fill_others(self, mol: Chem.Mol, distance_matrix: np.array, good_indices: List[int]):
        """
        Nan fill the inidices that are not the good_indices.

        :param mol:
        :param distance_matrix:
        :param good_indices:
        :return:
        """
        others = np.array(list(set(range(mol.GetNumAtoms())).difference(good_indices)))
        distance_matrix[others, :] = np.nan
        distance_matrix[:, others] = np.nan

    def _get_joining_penalties(self, combo: Chem.Mol, shape: Tuple[int, int]) -> np.ndarray:
        """
        Called by ``_find_closest``.
        THis is different from _get_merging_penalties

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

    def _nan_fill_submatrix(self, matrix: np.ndarray, indices: List[int]) -> None:
        """
        Given a square matrix, blank the self-submatrix of the group of indices
        There is probably a better way to do this.
        changed from _nan_submatrix as to nan is not a verb.

        :param matrix:
        :param indices:
        :return:
        """
        dimension = matrix.shape[0]
        bool_vector = np.zeros((dimension, 1)).astype(bool)
        indices = [i for i in indices if isinstance(i, (int, np.int64)) and i < dimension]
        bool_vector[indices] = True
        bool_matrix = np.tile(bool_vector, (1, dimension))
        logic = np.logical_and(bool_matrix, bool_matrix.transpose())
        matrix[logic] = np.nan

    # ============= Deletion ===========================================================================================

    def _mark_for_deletion(self, mol: Chem.Mol, i: int):
        mol.GetAtomWithIdx(i).SetBoolProp('DELETE', True)

    def _delete_marked(self, mol: Chem.RWMol):
        morituri = list(reversed(mol.GetAtomsMatchingQuery(Chem.rdqueries.HasPropQueryAtom('DELETE'))))
        for atom in morituri:
            mol.RemoveAtom(atom.GetIdx())

    # ============= Other ==============================================================================================

    def _copy_bonding(self, mol, keeper_idx: int, reject_idx: int, force: Optional[bool] = None):
        """
        formerly called `absorb`. Preps for absorbing. remove reject_idx separately.
        So copy bonding from reject_idx to keeper_idx.

        :param mol:
        :param keeper_idx:
        :param reject_idx:
        :param force:
        :return:
        """
        self.journal.debug(f'Copying atom bonding from {reject_idx} to {keeper_idx}, force={force}')
        keeper = mol.GetAtomWithIdx(keeper_idx)
        reject = mol.GetAtomWithIdx(reject_idx)
        # prevent confusion when it comes to triangles
        vertex = self._get_triangle(reject, keeper)
        if vertex:
            mol.RemoveBond(reject_idx, vertex)
        # deal with neighbours of reject atom
        for neighbor in reject.GetNeighbors():
            # collect bonding details between neighbour of reject and reject itself
            neigh_idx = neighbor.GetIdx()
            old_bond = mol.GetBondBetweenAtoms(reject_idx, neigh_idx)
            bt = old_bond.GetBondType()
            # forcing?
            if force is not None:
                force_bond = bool(force)
            elif old_bond.HasProp('_IsRingBond'):
                # the ring needs to be force!
                force_bond = True
            else:
                force_bond = False
            # mol.RemoveBond(j, n)
            if keeper_idx == neigh_idx:  # the neighbour is the keeper
                continue
            else:
                # copy bond. The provenance should be 'other_novel' not the original vai
                # provenance = BondProvenance.get_bond(old_bond)
                if force_bond and mol.GetBondBetweenAtoms(keeper_idx, neigh_idx) is None:
                    self.journal.debug(f'... Forcing bond between {keeper_idx} and {neigh_idx} during bond copying')
                    self._add_bond_regardlessly(mol=mol,
                                                first=keeper,
                                                second=neighbor,
                                                bond_type=bt,
                                                provenance='other_novel'
                                                )
                else:
                    self.journal.debug(
                        f'... Potentially adding bond between {keeper_idx} and {neigh_idx} during bond copying')
                    self._add_bond_if_possible(mol=mol,
                                               first=keeper,
                                               second=neighbor,
                                               provenance='other_novel')

    def _add_bond_regardlessly(self, mol, first: Chem.Atom, second: Chem.Atom, bond_type, provenance='other_novel'):
        """
        This methods does no checking and operates dangerously!

        :param mol:
        :param first:
        :param second:
        :param bond_type:
        :param provenance:
        :return:
        """
        first_idx = first.GetIdx()
        second_idx = second.GetIdx()
        # add if absent... (error prevention)
        present_bond = mol.GetBondBetweenAtoms(first_idx, second_idx)
        if present_bond is None:
            mol.AddBond(first_idx, second_idx, bond_type)
        new_bond = mol.GetBondBetweenAtoms(first_idx, second_idx)
        BondProvenance.set_bond(new_bond, provenance)

    def _assess_atom_for_possible_bonding(self, atom: Chem.Atom, bt: Chem.BondType) -> bool:
        """
        Method for add_bond_if_possible
        True means add, False means delete

        :param atom:
        :param bt:
        :return:
        """
        n_neigh = sum([self._is_count_valid(neigh) for neigh in atom.GetNeighbors()])
        if atom.GetAtomicNum() > 8:
            return True
        elif atom.HasProp('DELETE'):  # if it is to be deleted it should be fine.
            return True
        elif n_neigh <= 2 and atom.GetIsAromatic():
            return True  # Chem.BondType.SINGLE
        elif n_neigh <= 3 and not atom.GetIsAromatic():  # noqa Chem.Atom.GetIsAromatic(atom)
            return True
        else:
            return False  # too bonded already!

    def _add_bond_if_possible(self, mol, first: Chem.Atom, second: Chem.Atom, provenance='other_novel'):
        """
        This method is used by _copy_bonding, but triggered when force=False

        :param mol:
        :param first:
        :param second:
        :param provenance:
        :return:
        """
        # --- Prep
        first_idx = first.GetIdx()
        second_idx = second.GetIdx()
        present_bond = mol.GetBondBetweenAtoms(first_idx, second_idx)
        if second.HasProp('_ring_bond'):
            bt = getattr(Chem.BondType, second.GetProp('_ring_bond'))
        else:
            bt = None
        # === Bond already exists series ============================================
        # --- Bond already exists and not bond type is specified
        if present_bond is not None and bt is None:
            self.journal.debug(f'Bond between {first_idx} and {second_idx} already exists')
            return True  # exists
        # --- Bond already exists but has an error
        elif present_bond is not None and present_bond.GetBondType() is None:
            present_bond.SetBondType(Chem.BondType.SINGLE)
            return True
        # --- Bond already exists and matches the expected bond type
        elif present_bond is not None and bt.name == present_bond.GetBondType().name:
            return True
        # --- Bond already exists but does not match the expected bond type
        elif present_bond is not None:
            present_bond.SetBondType(bt)
            return True
        # === Assess if new bond should be added ============================================
        # --- Don't add if it would make a triangle
        elif self._is_would_be_triangle(first, second):
            self.journal.debug(f'Bond between {first_idx} and {second_idx} would make a triangle, skipping')
            return False
        # --- Don't add if it would make a square
        elif self._is_would_be_square(first, second):
            self.journal.debug(f'Bond between {first_idx} and {second_idx} would make a square, skipping')
            return False
        # --- Don't add if it would ruin the warhead
        elif self._is_connected_warhead(second, first):
            self.journal.debug(f'Bond between {first_idx} and {second_idx} would break a warhead, skipping')
            return False
        # --- New bond is green lit
        elif self._assess_atom_for_possible_bonding(first, bt) and self._assess_atom_for_possible_bonding(second, bt):
            mol.AddBond(first_idx, second_idx, bt if bt is not None else Chem.BondType.SINGLE)
            new_bond = mol.GetBondBetweenAtoms(first_idx, second_idx)
            BondProvenance.set_bond(new_bond, provenance)
            return True
        # --- New bond is no go
        else:
            # len(Chem.GetMolFrags(mol, sanitizeFrags=False)) ought to be checked.
            # however, join_neighboring gets called by emergency bonding so should be fine.
            return False  # too bonded etc.

    # === conditional selectors ========================================================================================

    def _is_would_be_triangle(self, first: Chem.Atom, second: Chem.Atom) -> bool:
        """
        Get bool of whether two atoms share a common neighbor. Ie. joining them would make a triangle.
        Direct bond does not count.

        :param first:
        :param second:
        :return:
        """
        if self._get_triangle(first, second) is not None:
            return True
        else:
            return False

    def _is_would_be_square(self, first: Chem.Atom, second: Chem.Atom) -> bool:
        """
        Get bool of whether two atoms share a common neighbor+over-neighbor. Ie. joining them would make a square.
        Direct bond does not count.

        :param first:
        :param second:
        :return:
        """
        for third in [neigh for neigh in second.GetNeighbors() if neigh.GetIdx() != first.GetIdx()]:
            if self._is_would_be_triangle(first, third) is True:
                return True
        else:
            return False

    def _get_triangle(self, first: Chem.Atom, second: Chem.Atom) -> Union[int, None]:
        """
        Get the third atom...

        :param first: atom
        :param second: atom
        :return: atom index of third
        """
        triang = self._get_triangles(first, second)
        if triang:
            return triang[0]
        else:
            return None

    def _get_triangles(self, first: Chem.Atom, second: Chem.Atom) -> Union[List[int], None]:
        """
        Get the third atoms... (Square situation)

        :param first: atom
        :param second: atom
        :return: atom index of third
        """
        get_neigh_idxs = lambda atom: [neigh.GetIdx() for neigh in atom.GetNeighbors() if
                                       self._is_count_valid(neigh)]
        f_neighs: List[int] = get_neigh_idxs(first)
        s_neighs: List[int] = get_neigh_idxs(second)
        a = set(f_neighs) - {first.GetIdx(), second.GetIdx()}
        b = set(s_neighs) - {first.GetIdx(), second.GetIdx()}
        others: List[int] = list(a.intersection(b))
        if len(others) == 0:  # is a disjoined
            return None
        else:
            return list(others)

    def _is_count_valid(self, atom: Chem.Atom) -> bool:
        """
        Some atoms are not to be counted as they will be deleted.

        :param atom:
        :return:
        """
        if atom.HasProp('DELETE'):
            return False
        elif atom.HasProp('_ori_i') and atom.GetIntProp('_ori_i') == -1:
            return False
        else:
            return True

    def _get_square(self, first: Chem.Atom, second: Chem.Atom) -> Union[Tuple[int, int], None]:
        for third in [neigh for neigh in second.GetNeighbors() if neigh.GetIdx() != first.GetIdx()]:
            fourths = self._get_triangles(first, third)
            if fourths and len(fourths) > 1:
                fourth = [f for f in fourths if f != second.GetIdx()][0]
                return third.GetIdx(), fourth
        else:
            return None

    def _is_connected_warhead(self, atom, anchor_atom):
        if not atom.HasProp('_Warhead'):
            return False
        elif atom.GetBoolProp('_Warhead') == False:
            return False
        else:
            # is it a single compound?
            frags = Chem.GetMolFrags(atom.GetOwningMol(), sanitizeFrags=False)
            if len(frags) == 1:
                return True
            else:
                for frag in frags:
                    if atom.GetIdx() in frag and anchor_atom.GetIdx() in frag:
                        return True
                    elif atom.GetIdx() in frag:
                        return False  # if the warhead is not connected pretend it is not a warhead.
                    else:
                        pass
                else:
                    raise ValueError('I do not think this is possible.')
