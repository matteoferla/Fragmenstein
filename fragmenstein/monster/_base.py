import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign, rdmolops
from rdkit.Geometry.rdGeometry import Point3D
from typing import Optional, Dict, List, Any, Tuple, Union
from .bond_provenance import BondProvenance
import logging

class _MonsterBaseMixin:

    journal = logging.getLogger('Fragmenstein')

    dummy_symbol = '*'
    dummy = Chem.MolFromSmiles(dummy_symbol)  #: The virtual atom where the targets attaches
    cutoff = 2.
    joining_cutoff = 5.  # how distant (in Å) is too much?
    atoms_in_bridge_cutoff = 2  # how many bridge atoms can be deleted? (0 = preserves norbornane, 1 = preserves monsterantane)
    throw_on_discard = False
    matching_modes = [
        dict(atomCompare=rdFMCS.AtomCompare.CompareAny,
             bondCompare=rdFMCS.BondCompare.CompareAny,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=False),  # this shape based matching is too permissive,
        dict(atomCompare=rdFMCS.AtomCompare.CompareAny,
             bondCompare=rdFMCS.BondCompare.CompareOrder,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=False),
        dict(atomCompare=rdFMCS.AtomCompare.CompareElements,
             bondCompare=rdFMCS.BondCompare.CompareOrder,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=False),
        dict(atomCompare=rdFMCS.AtomCompare.CompareAny,
             bondCompare=rdFMCS.BondCompare.CompareAny,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=True),
        dict(atomCompare=rdFMCS.AtomCompare.CompareAny,
             bondCompare=rdFMCS.BondCompare.CompareOrder,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=True),
        dict(atomCompare=rdFMCS.AtomCompare.CompareElements,
             bondCompare=rdFMCS.BondCompare.CompareOrder,
             ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
             ringMatchesRingOnly=True)]

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
        return atom.HasProp('_Warhead') and atom.GetBoolProp('_Warhead') is True

    # func: https://stackoverflow.com/questions/41921255/staticmethod-object-is-not-callable
    closeness_weights =[
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
        penalties = self._get_penalties(combo, distance_matrix.shape)
        # ========= get closest
        pendist_matrix = penalties + distance_matrix
        pendistance = np.nanmin(pendist_matrix)
        if np.isnan(pendistance):
            raise ConnectionError('This is impossible. Previous is absent??')
        else:
            candidates = []

            def get_closest(pendistance):
                p = np.where(pendist_matrix == pendistance)
                anchor_A = int(p[0][0])
                anchor_B = int(p[1][0])
                distance = distance_matrix[anchor_A, anchor_B]
                penalty = penalties[anchor_A, anchor_B]
                self.journal.debug(f'Connecting {anchor_A} with {anchor_B}, {penalty} penalised distance of {distance}')
                return anchor_A, anchor_B, distance

            anchor_A, anchor_B, distance = get_closest(pendistance)
            candidates.append((anchor_A, anchor_B, distance))
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


    def _get_distance_matrix(self, combo: Chem.Mol, A: Union[Chem.Mol, np.ndarray], B: Union[Chem.Mol, np.ndarray]) -> np.ndarray:
        """
        Called by ``_find_closest`` and ``_determine_mergers_novel_ringcore_pair`` in collapse ring (for expansion).

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

    # ============= Deletion ===========================================================================================

    def _mark_for_deletion(self, mol: Chem.Mol, i: int):
        mol.GetAtomWithIdx(i).SetBoolProp('DELETE', True)

    def _delete_marked(self, mol: Chem.RWMol):
        morituri = list(reversed(mol.GetAtomsMatchingQuery(Chem.rdqueries.HasPropQueryAtom('DELETE'))))
        for atom in morituri:
            mol.RemoveAtom(atom.GetIdx())

    # ============= Other ==============================================================================================

    def _copy_bonding(self, mol, i: int, j: int, force: Optional[bool] = None):
        """
        formerly called `absorb`. Preps for absorbing. remove J separately.

        :param mol:
        :param i:
        :param j:
        :return:
        """
        self.journal.debug(f'Absorbing atom {i} with {j}')
        absorbenda = mol.GetAtomWithIdx(i)
        absorbiturum = mol.GetAtomWithIdx(j)
        for neighbor in absorbiturum.GetNeighbors():
            neigh_i = neighbor.GetIdx()
            old_bond = mol.GetBondBetweenAtoms(j, neigh_i)
            bt = old_bond.GetBondType()
            if force is not None:
                pass
            if old_bond.HasProp('_IsRingBond'):
                # the ring needs to be force!
                force = True
            else:
                force = False
            # mol.RemoveBond(j, n)
            if i == neigh_i:
                continue
            else:
                atom_i, atom_j = mol.GetAtomWithIdx(i), mol.GetAtomWithIdx(j)
                if force and mol.GetBondBetweenAtoms(i, neigh_i) is None:
                    self.journal.debug(f'Forcing bond between {i} and {neigh_i}')
                    mol.AddBond(i, neigh_i, bt)
                    new_bond = mol.GetBondBetweenAtoms(i, neigh_i)
                    BondProvenance.copy_bond(old_bond, new_bond)
                else:
                    self._add_bond_if_possible(mol, atom_i, atom_j)

    def _add_bond_if_possible(self, mol, atom_i, atom_j, provenance='other_novel'):
        i = atom_i.GetIdx()
        j = atom_j.GetIdx()

        def assess_atom(atom: Chem.Atom, bt: Chem.BondType) -> Tuple[bool, Chem.BondType]:
            """
            True means add, False means delete

            :param atom:
            :param bt:
            :return:
            """
            n_neigh = sum([self._is_count_valid(neigh) for neigh in atom.GetNeighbors()])
            if atom.GetAtomicNum() > 8:
                return True, bt
            elif atom.HasProp('DELETE'):  # if it is to be deleted it should be fine.
                return True, bt
            elif n_neigh <= 2 and atom.GetIsAromatic():
                return True, Chem.BondType.SINGLE
            elif n_neigh <= 3 and not atom.GetIsAromatic():
                return True, bt
            else:
                return False, bt  # too bonded already!

        if self._is_triangle(atom_i, atom_j):
            self.journal.debug(f'Bond between {i} and {j} would make a triangle, skipping')
            return False
        elif self._is_square(atom_i, atom_j):
            self.journal.debug(f'Bond between {i} and {j} would make a square, skipping')
            return False
        elif self._is_connected_warhead(atom_j, atom_i):
            self.journal.debug(f'Bond between {i} and {j} would break a warhead, skipping')
            return False
        else:
            present_bond = mol.GetBondBetweenAtoms(i, j)
            if atom_j.HasProp('_ring_bond'):
                bt = getattr(Chem.BondType, atom_j.GetProp('_ring_bond'))
            else:
                bt = None
            if present_bond is not None and bt is None:
                self.journal.debug(f'Bond between {i} and {j} already exists')
                pass  # exists
            elif present_bond is not None and present_bond.GetBondType() is None:
                present_bond.SetBondType(Chem.BondType.SINGLE)
            elif present_bond is not None and present_bond.GetBondType() is not None and bt.name == present_bond.GetBondType().name:
                pass  # exists and has correct bond
            elif present_bond is not None:
                present_bond.SetBondType(bt)
                return True
            else:
                v, bt = assess_atom(atom_i, bt)
                w, bt = assess_atom(atom_j, bt)
                if v and w and bt is not None:
                    mol.AddBond(i, j, bt)
                    new_bond = mol.GetBondBetweenAtoms(i, j)
                    BondProvenance.set_bond(new_bond, provenance)
                    return True
                elif v and w:
                    mol.AddBond(i, j, Chem.BondType.SINGLE)
                    new_bond = mol.GetBondBetweenAtoms(i, j)
                    BondProvenance.set_bond(new_bond, provenance)
                    return True
                else:
                    # len(Chem.GetMolFrags(mol, sanitizeFrags=False)) ought to be checked.
                    # however, join_neighboring gets called by emergency bonding so should be fine.
                    return False  # too bonded etc.

    # === conditional selectors ========================================================================================

    def _is_triangle(self, first: Chem.Atom, second: Chem.Atom) -> bool:
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
        f_neighs = get_neigh_idxs(first)
        s_neighs = get_neigh_idxs(second)
        a = set(f_neighs) - {first.GetIdx(), second.GetIdx()}
        b = set(s_neighs) - {first.GetIdx(), second.GetIdx()}
        others = list(a.intersection(b))
        if len(others) == 0:  # is disjoined
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

    def _is_square(self, first: Chem.Atom, second: Chem.Atom) -> bool:
        """
        Get bool of whether two atoms share a common over-neighbor. Ie. joining them would make a square.
        Direct bond does not count.

        :param first:
        :param second:
        :return:
        """
        for third in [neigh for neigh in second.GetNeighbors() if neigh.GetIdx() != first.GetIdx()]:
            if self._is_triangle(first, third) is True:
                return True
        else:
            return False

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
