########################################################################################################################

__doc__ = \
    """
This is the ring collapsing code.
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################

import json, itertools
from warnings import warn
from rdkit.Geometry.rdGeometry import Point3D
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Optional, Dict, List, Any, Tuple, Union, Callable
import numpy as np
from collections import Counter
from functools import partial
from .bond_provenance import BondProvenance
from ._base import _AdamBaseMixin


class _AdamRing(_AdamBaseMixin):
    def __init__(self, _debug_draw=False):
        # abstracted...
        self._debug_draw = _debug_draw

    def collapse_mols(self, mols: List[Chem.Mol]):
        mols = [self.collapse_ring(mol) for mol in mols]
        [self.offset(mol) for mol in mols]
        return mols

    def store_positions(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Saves positional data as _x, _y, _z and majorly ``_ori_i``, the original index.
        The latter gets used by ``_get_new_index``.

        :param mol:
        :return:
        """
        conf = mol.GetConformer()
        name = mol.GetProp('_Name')
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atom.SetIntProp('_ori_i', i)
            atom.SetProp('_ori_name', name)
            atom.SetDoubleProp('_x', pos.x)
            atom.SetDoubleProp('_y', pos.y)
            atom.SetDoubleProp('_z', pos.z)
        return mol


    # =========== Collapse & Expand ====================================================================================

    def collapse_ring(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Collapses a ring(s) into a single dummy atom(s).
        Stores data as JSON in the atom.

        :param mol:
        :return:
        """
        self.store_positions(mol)
        mol = Chem.RWMol(mol)
        conf = mol.GetConformer()
        center_idxs = []
        morituri = []
        old2center = defaultdict(list)
        # Store simple info
        for atomset in mol.GetRingInfo().AtomRings():
            morituri.extend(atomset)
            neighs = []
            neighbonds = []
            bonds = []
            xs = []
            ys = []
            zs = []
            elements = []
            # add elemental ring
            c = mol.AddAtom(Chem.Atom('C'))
            center_idxs.append(c)
            central = mol.GetAtomWithIdx(c)
            name = mol.GetProp('_Name') if mol.HasProp('_Name') else '???'
            central.SetProp('_ori_name', name),
            # get data for storage
            for i in atomset:
                old2center[i].append(c)
                atom = mol.GetAtomWithIdx(i)
                neigh_i = [a.GetIdx() for a in atom.GetNeighbors()]
                neighs.append(neigh_i)
                bond = [mol.GetBondBetweenAtoms(i, j).GetBondType().name for j in neigh_i]
                bonds.append(bond)
                pos = conf.GetAtomPosition(i)
                xs.append(pos.x)
                ys.append(pos.y)
                zs.append(pos.z)
                elements.append(atom.GetSymbol())
            # store data in elemental ring
            central.SetIntProp('_ori_i', -1)
            central.SetProp('_ori_is', json.dumps(atomset))
            central.SetProp('_neighbors', json.dumps(neighs))
            central.SetProp('_xs', json.dumps(xs))
            central.SetProp('_ys', json.dumps(ys))
            central.SetProp('_zs', json.dumps(zs))
            central.SetProp('_elements', json.dumps(elements))
            central.SetProp('_bonds', json.dumps(bonds))
            conf.SetAtomPosition(c, Point3D(*[sum(axis) / len(axis) for axis in (xs, ys, zs)]))
        # Store complex info
        for atomset, center_i in zip(mol.GetRingInfo().AtomRings(), center_idxs):
            # bond to elemental ring
            central = mol.GetAtomWithIdx(center_i)
            neighss = json.loads(central.GetProp('_neighbors'))
            bondss = json.loads(central.GetProp('_bonds'))
            for neighs, bonds in zip(neighss, bondss):
                for neigh, bond in zip(neighs, bonds):
                    if neigh not in atomset:
                        bt = getattr(Chem.BondType, bond)
                        if neigh not in morituri:
                            mol.AddBond(center_i, neigh, bt)
                            new_bond = mol.GetBondBetweenAtoms(center_i, neigh)
                            BondProvenance.set_bond(new_bond, 'original')
                        else:
                            for other_center_i in old2center[neigh]:
                                if center_i != other_center_i:
                                    if not mol.GetBondBetweenAtoms(center_i, other_center_i):
                                        mol.AddBond(center_i, other_center_i, bt)
                                        new_bond = mol.GetBondBetweenAtoms(center_i, other_center_i)
                                        BondProvenance.set_bond(new_bond, 'original')
                                    break
                            else:
                                raise ValueError(f'Cannot find what {neigh} became')
        # Remove atoms
        for i in sorted(set(morituri), reverse=True):
            mol.RemoveAtom(self._get_new_index(mol, i))
        return mol.GetMol()

    def expand_ring(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Undoes collapse ring

        :param mol: untouched.
        :return:
        """
        self.journal.debug('Starting ring expansion')
        mol = Chem.RWMol(mol)
        rings = self._get_expansion_data(mol)  # List[Dict[str, List[Any]]]
        self._place_ring_atoms(mol, rings)
        # bonded_as_original. Rectifier will fix.
        self._restore_original_bonding(mol, rings)
        self._add_novel_bonding(mol, rings) # formerly `_ring_overlap_scenario` and `_infer_bonding_by_proximity`.
        self._delete_collapsed(mol)
        self._detriangulate(mol)
        mol = self._emergency_joining(mol) # does not modify in place!
        if mol is None:
            raise ValueError('(Impossible) Failed at some point...')
        elif isinstance(mol, Chem.RWMol):
            return mol.GetMol()
        else:
            return mol

    # =========== Offset ===============================================================================================

    _collapsed_ring_offset = 0

    def offset(self, mol: Chem.Mol):
        """
        This is to prevent clashes.
        The numbers of the ori indices stored in collapsed rings are offset by the class variable (_collapsed_ring_offset)
        multiples of 100. (autoincrements to avoid dramas)

        :param mol:
        :return:
        """
        self._collapsed_ring_offset += 100
        old2new = {}
        # sort the not ringcore atoms
        for atom in mol.GetAtoms():
            if atom.GetIntProp('_ori_i') != -1:
                o = atom.GetIntProp('_ori_i')
                n = o + self._collapsed_ring_offset
                atom.SetIntProp('_ori_i', n)
                old2new[o] = n
            else:
                pass # ringcores have -1 ori_i
        # sort the ringcore
        for atom in self._get_collapsed_atoms(mol):
            old = json.loads(atom.GetProp('_ori_is'))
            new = [i + self._collapsed_ring_offset for i in old]
            atom.SetProp('_ori_is', json.dumps(new))
            old2new = {**old2new, **dict(zip(old, new))}
        # this has to be done afterwards in case of a bonded mol
        for atom in self._get_collapsed_atoms(mol):
            old_neighss = json.loads(atom.GetProp('_neighbors')) #if i in old2new else i
            new_neighss = [[old2new[i] for i in old_neighs if i in old2new] for old_neighs in old_neighss]
            atom.SetProp('_neighbors', json.dumps(new_neighss))

    def _renumber_original_indices(self, mol: Chem.Mol,
                                   mapping: Dict[int, int],
                                   name_restriction: Optional[str] = None):
        """
        Renumbers the indices in the ring virtual atom to new mapping.
        this is unused because the code needing it unfinished.

        :param mol:
        :param mapping: old index to new index dictionary
        :param name_restriction:
        :return:
        """
        for atom in mol.GetAtoms():
            if name_restriction is not None and \
                    atom.HasProp('_ori_name') and \
                    atom.GetProp('_ori_name') != name_restriction:
                continue
            i = atom.GetIntProp('_ori_i')
            if i == -1:
                # original indices
                ori = json.loads(atom.GetProp('_ori_is'))
                alt = [dd if dd not in mapping else mapping[dd] for dd in ori]
                atom.SetProp('_ori_is', json.dumps(alt))
                # original neighbors
                ori = json.loads(atom.GetProp('_neighbors'))
                alt = [[dd if dd not in mapping else mapping[dd] for dd in inner] for inner in ori]
                atom.SetProp('_neighbors', json.dumps(alt))
            elif i in mapping:
                atom.SetIntProp('_ori_i', mapping[i])
            else:
                pass

    # =========== Expand data ==========================================================================================

    def _get_expansion_data(self, mol: Chem.Mol) -> List[Dict[str, List[Any]]]:
        """
        Returns a list for each collapsed ring marking atom each with a dictionary.
        Example:

             {'atom': <rdkit.Chem.rdchem.Atom at 0x7f926fafb030>,
              'ori_name': 'APR',
              'elements': ['O', 'C', 'C', 'C', 'C'],
              'neighbors': [[4, 10], [2, 5, 6], [4, 7, 8], [6, 9, 10], [5, 8, 35]],
              'ori_is': [5, 4, 6, 8, 10],
              'xs': [0.55, 1.16, 0.785, -0.527, -0.198],
              'ys': [15.473, 15.587, 14.293, 13.909, 14.28],
              'zs': [-3.205, -4.516, -5.259, -4.577, -3.132],
              'bonds': [['SINGLE', 'SINGLE'],
               ['SINGLE', 'SINGLE', 'SINGLE'],
               ['SINGLE', 'SINGLE', 'SINGLE'],
               ['SINGLE', 'SINGLE', 'SINGLE'],
               ['SINGLE', 'SINGLE', 'SINGLE']]}

        :param mol:
        :return:
        """
        return [
            dict(atom=atom,
                 ori_name=atom.GetProp('_ori_name'),
                 elements=json.loads(atom.GetProp('_elements')),
                 neighbors=json.loads(atom.GetProp('_neighbors')),
                 ori_is=json.loads(atom.GetProp('_ori_is')),
                 xs=json.loads(atom.GetProp('_xs')),
                 ys=json.loads(atom.GetProp('_ys')),
                 zs=json.loads(atom.GetProp('_zs')),
                 bonds=json.loads(atom.GetProp('_bonds')))
            for atom in self._get_collapsed_atoms(mol)]

    def _get_expansion_for_atom(self, data: Dict[str, List[Any]], i: int) -> Dict[str, Any]:
        return {k.replace('s', ''): data[k][i] if isinstance(data[k], list) else data[k] for k in data}

    # === Key steps ====================================================================================================

    def _place_ring_atoms(self, mol, rings: List[Dict[str, Union[Chem.Atom, List[Any]]]]):
        """
        Plase all atoms stored in rings.

        :param mol:
        :param rings:
        :return:
        """
        conf = mol.GetConformer()
        for ring in rings:  # atoms in ring addition
            indices = [] # will store current indices
            for i in range(len(ring['elements'])):  # atom addition
                collapsed_atom_data = self._get_expansion_for_atom(ring, i)
                # atom will be added if it is not already present!
                if self._is_present(mol, collapsed_atom_data['ori_i']):
                    natom = self._get_new_index(mol, collapsed_atom_data['ori_i'], search_collapsed=False)
                    self.journal.debug(f"{natom} (formerly {collapsed_atom_data['ori_i']} existed already." +
                                  "Fused ring or similar.")
                else:
                    n = mol.AddAtom(Chem.Atom(collapsed_atom_data['element']))
                    natom = mol.GetAtomWithIdx(n)
                    conf.SetAtomPosition(n, Point3D(collapsed_atom_data['x'],
                                                    collapsed_atom_data['y'],
                                                    collapsed_atom_data['z']))
                    natom.SetIntProp('_ori_i', collapsed_atom_data['ori_i'])
                    natom.SetDoubleProp('_x', collapsed_atom_data['x'])
                    natom.SetDoubleProp('_y', collapsed_atom_data['y'])
                    natom.SetDoubleProp('_z', collapsed_atom_data['z'])
                    natom.SetProp('_ori_name', collapsed_atom_data['ori_name'])
                    indices.append(n)
            ringcore = ring['atom']
            ringcore.SetProp('_current_is', json.dumps(indices))

    def _restore_original_bonding(self, mol: Chem.RWMol, rings: List[Dict[str, List[Any]]]) -> None:
        self.journal.debug('Restoring original bonding if any.')
        to_be_waited_for = []
        for ring in rings:
            for i in range(len(ring['elements'])):
                # iteration per atom:
                collapsed_atom_data = self._get_expansion_for_atom(ring, i)
                old_i = collapsed_atom_data['ori_i']
                new_i = self._get_new_index(mol, old_i, search_collapsed=False)
                for old_neigh, bond in zip(collapsed_atom_data['neighbor'], collapsed_atom_data['bond']):
                    bt = getattr(Chem.BondType, bond)
                    new_neigh = self._get_new_index(mol, old_neigh, search_collapsed=False)
                    info = f'new_i={new_i}, new_neig={new_neigh}, old_i={old_i}, old_neigh={old_neigh}'
                    present_bond = mol.GetBondBetweenAtoms(new_i, new_neigh)
                    if present_bond is None:
                        assert new_i != new_neigh, f'Cannot bond to self. {info}'
                        mol.AddBond(new_i, new_neigh, bt)
                        present_bond = mol.GetBondBetweenAtoms(new_i, new_neigh)
                        present_bond.SetBoolProp('_IsRingBond', True)
                        BondProvenance.set_bond(present_bond, 'original')
                        distance = Chem.rdMolTransforms.GetBondLength(mol.GetConformer(), new_i, new_neigh)
                        assert distance < 4, f'Bond length too long ({distance}. {info}'
                    elif present_bond.GetBondType().name != bond:
                        self.journal.warning(f'bond between {new_i} {new_neigh} exists already ' +
                                    f'(has {present_bond.GetBondType().name} expected {bt})')
                        present_bond.SetBondType(bt)
                        present_bond.SetBoolProp('_IsRingBond', True)
                        BondProvenance.set_bond(present_bond, 'original')
                    else:
                        self.journal.debug(f'bond between {new_i} {new_neigh} exists already ' +
                                  f'(has {present_bond.GetBondType().name} expected {bt})')
                        pass
        #             try:
        #
        #             except ValueError:
        #                 log.warning(f"The neighbour {old_neigh} of {collapsed_atom_data['ori_i']} with {bt} " +
        #                             "does not yet exist")
        #                 to_be_waited_for.append((new_i, old_neigh, bt))
        # for new_i, old_neigh, bt in to_be_waited_for:
        #     try:
        #         new_neigh = self._get_new_index(mol, old_neigh,
        #                                         name_restriction=mol.GetAtomWithIdx(new_i).GetProp('_ori_name'))
        #         print(f'{old_neigh} was missing, but has appeared since as {new_neigh}')
        #         if not mol.GetBondBetweenAtoms(new_i, new_neigh):
        #             add_bond(mol, new_i, new_neigh, bt)
        #     except (KeyError, ValueError) as err:
        #         warn(str(err))

    def _delete_collapsed(self, mol: Chem.RWMol):
        for a in reversed(range(mol.GetNumAtoms())):
            if mol.GetAtomWithIdx(a).GetIntProp('_ori_i') == -1:
                mol.RemoveAtom(a)

    def _add_novel_bonding(self, mol: Chem.RWMol, rings: List[Dict[str, List[Any]]]):
        """
        Formerly `_ring_overlap_scenario`,`_connenct_ring`,  `_infer_bonding_by_proximity`.

        Two close rings can be one of:

        * bonded rings: Two separate rings bonded
        * Spiro rings: One atom in common
        * fused rings: Two atoms in common
        * bridged rings: counts as two rings, not three.

        The atom placement stops the same atom being placed twice...
        ...But the major issue is that the rings may and often do come from two different hits.
        Hence why it does not store it in memory.

        :param mol:
        :param rings: output of `_get_expansion_data`.
        :type rings: List[Dict[str, List[Any]]]
        :return:
        """
        # ===== Deal with Ring on ring bonding
        self.journal.debug('Adding novel bonding (if any)...')
        novel_ringcore_pairs = self._get_novel_ringcore_pairs(rings)
        for ringcore_A, ringcore_B in novel_ringcore_pairs:
            self.journal.debug(f'determining novel bond between {ringcore_A} and {ringcore_B}')
            # _determine_mergers_novel_ringcore_pair finds mergers
            self._determine_mergers_novel_ringcore_pair(mol, ringcore_A, ringcore_B)
        # ===== Deal with Ring on other bonding
        # formerly: _infer_bonding_by_proximity
        novel_other_pairs = self._get_novel_other_pairs(rings)
        for ringcore, other in novel_other_pairs:
            self.journal.debug(f'determining novel bond between {ringcore} and {other}')
            # _determine_mergers_novel_ringcore_pair finds, bonds and marks for deletion.
            self._determine_mergers_novel_other_pair(mol, ringcore, other)
        # ===== Clean up
        self._delete_marked(mol)

    # =========== dependant methods =================================================================================

    def _get_new_index(self, mol: Chem.Mol, old: int, search_collapsed=True,
                       name_restriction: Optional[str] = None) -> int:
        """
        Given an old index check in ``_ori_i`` for what the current one is.
        NB. ring placeholder will be -1 and these also have ``_ori_is``. a JSON of the ori_i they summarise.

        :param mol:
        :param old: old index
        :param search_collapsed: seach also in ``_ori_is``
        :parm name_restriction: restrict to original name.
        :return:
        """
        for i, atom in enumerate(mol.GetAtoms()):
            if name_restriction is not None and atom.GetProp('_ori_name') != name_restriction:
                pass
            elif atom.GetIntProp('_ori_i') == old:
                return i
            elif search_collapsed and \
                    atom.HasProp('_ori_is') and \
                    old in json.loads(atom.GetProp('_ori_is')):
                return i
            else:
                pass
        else:
            raise ValueError(f'Cannot find {old}')

    def _is_present(self, mol, i):
        """
        Find if in ``mol`` there is an atom whose original index was ``i``.
        Wrapper around ``_get_new_index``, but does not return the index.

        :param mol:
        :param i: original index
        :return:
        """
        try:
            self._get_new_index(mol, i, search_collapsed=False)
            # raises value error.
            return True
        except ValueError as err:  # no atom is present (actually the default)
            return False

    def _get_collapsed_atoms(self, mol: Chem.Mol) -> List[Chem.Atom]:
        return [atom for atom in mol.GetAtoms() if atom.GetIntProp('_ori_i') == -1]

    # ==== Novel Ring to Ring bonding ==================================================================================

    def _get_novel_ringcore_pairs(self, rings) -> List[Tuple[Chem.Atom, Chem.Atom]]:
        pairs = []
        for ring in rings:
            # ring: Dict[str, List[Any]]
            ringcore = ring['atom']  # Chem.Atom
            for neigh in ringcore.GetNeighbors():
                # find those ringcore atoms that are connected. via these conditions:
                has_ringcore_neighbor = neigh.HasProp('_ori_name') and neigh.GetIntProp('_ori_i') == -1
                is_new_pair = ringcore.GetIdx() < neigh.GetIdx()  # to avoid doing it twice each direction
                is_novel_connection = neigh.HasProp('_ori_name') and ring['ori_name'] != neigh.GetProp('_ori_name')
                if has_ringcore_neighbor and is_new_pair and is_novel_connection:
                    # This ringcore atom shares a novel border with another ringcore atom
                    pairs.append((ringcore, neigh))
        return pairs

    def _determine_mergers_novel_ringcore_pair(self,
                                               mol: Chem.RWMol,
                                               ringcore_A: Chem.Atom,
                                               ringcore_B: Chem.Atom) -> List[Tuple[int, int]]:
        """
        Formerly part of ``_ring_overlap_scenario``.
        Preps without deleting. ``_delete_marked(mol)`` does that.
        bonded, spiro, fused

        :param mol:
        :param ringcore_A:
        :param ringcore_B:
        :return: list of atoms to be merged
        """
        absorption_distance = 1. # Å
        indices_A = json.loads(ringcore_A.GetProp('_current_is'))
        indices_B = json.loads(ringcore_B.GetProp('_current_is'))
        distance_matrix = self._get_distance_matrix(mol, indices_A, indices_B) # currently in `_join_neighboring`.
        # distance matrix is for the whole thing
        # TODO merge into _get_distance_matrix
        # getting the other atom indices:
        others = np.array(list(set(range(mol.GetNumAtoms())).difference(indices_A + indices_B)))
        distance_matrix[others, :] = np.nan
        distance_matrix[:, others] = np.nan
        # get closest pair.
        distance = np.nanmin(distance_matrix)
        if np.isnan(distance):
            self.journal.critical('This is impossible. Two neighbouring rings cannot be connected.')
            return []
        elif distance > absorption_distance: # bonded
            p = np.where(distance_matrix == distance)
            a = int(p[0][0])
            b = int(p[1][0])
            present_bond = mol.GetBondBetweenAtoms(ringcore_A.GetIdx(), ringcore_B.GetIdx())
            bt = present_bond.GetBondType()
            if bt is None or bt == Chem.BondType.UNSPECIFIED:
                bt = Chem.BondType.SINGLE
            n = mol.AddBond(a, b, bt)
            new_bond = mol.GetBondBetweenAtoms(a, b)
            BondProvenance.copy_bond(present_bond, new_bond)
            self.journal.info('A novel bond-connected ring pair was found')
            self._mark_for_deletion(mol, b)
            self._copy_bonding(mol, a, b, force=True)
            return [] #bonded
        else:
            p = np.where(distance_matrix == distance)
            a = int(p[0][0])
            b = int(p[1][0])
            self._mark_for_deletion(mol, b)
            self._copy_bonding(mol, a, b, force=True)
            distance_matrix[a, b] = np.nan
            distance_matrix[b, a] = np.nan
            distance = np.nanmin(distance_matrix)
            if np.isnan(distance) or distance < absorption_distance:
                p = np.where(distance_matrix == distance)
                c = int(p[0][0])
                d = int(p[1][0])
                self._mark_for_deletion(mol, d)
                self._copy_bonding(mol, c, d, force=True)
                self.journal.info('A novel fused ring pair was found')
                return [(a, b), (c, d)] #fused
            else:
                self.journal.info('A novel spiro ring pair was found')
                return [(a,b)] # spiro

    # ==== Novel Ring to Ring bonding ==================================================================================

    # def _infer_bonding_by_proximity(self, mol):
    #     raise Exception
    #     # fix neighbours
    #     # this should not happen. But just in case!
    #     while True:
    #         ringsatoms = self._get_ring_info(mol)
    #         for ringA, ringB in itertools.combinations(ringsatoms, 2):
    #             n = mol.GetNumAtoms()
    #             self.absorb_overclose(mol, ringA, ringB, cutoff=1.)
    #             if n != mol.GetNumAtoms():
    #                 break
    #         else:
    #             break
    #     new_ringers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp('expanded')]
    #     self.absorb_overclose(mol, new_ringers)
    #     new_ringers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp('expanded')]
    #     self.join_overclose(mol, new_ringers)
    #     self.join_rings(mol)
    #     self._triangle_warn(mol)


    def _get_novel_other_pairs(self, rings) -> List[Tuple[Chem.Atom, Chem.Atom]]:
        ## similar to novel ringcore..but opposite
        pairs = []
        for ring in rings:
            # ring: Dict[str, List[Any]]
            ringcore = ring['atom']  # Chem.Atom
            for neigh in ringcore.GetNeighbors():
                # find those ringcore-other atoms that are connected. via these conditions:
                has_notringcore_neighbor = (not neigh.HasProp('_ori_name')) or neigh.GetIntProp('_ori_i') != -1
                is_novel_connection = (not neigh.HasProp('_ori_name')) or ring['ori_name'] != neigh.GetProp('_ori_name')
                if has_notringcore_neighbor and is_novel_connection:
                    # This ringcore atom shares a novel border with another ringcore atom
                    pairs.append((ringcore, neigh))
        return pairs

    def _determine_mergers_novel_other_pair(self,
                                               mol: Chem.RWMol,
                                               ringcore: Chem.Atom,
                                               other: Chem.Atom) -> List[Tuple[int, int]]:
        absorption_distance = 1.  # Å
        indices_ring = json.loads(ringcore.GetProp('_current_is'))
        indices_other = [other.GetIdx()]
        distance_matrix = self._get_distance_matrix(mol, indices_ring, indices_other)  # currently in `_join_neighboring`.
        # distance matrix is for the whole thing
        # TODO merge into _get_distance_matrix
        # getting the other atom indices:
        others = np.array(list(set(range(mol.GetNumAtoms())).difference(indices_ring + indices_other)))
        distance_matrix[others, :] = np.nan
        distance_matrix[:, others] = np.nan
        # penalties
        penalties = np.zeros(distance_matrix.shape)
        for i in indices_ring:
            atom = mol.GetAtomWithIdx(i)
            neighs = [neigh for neigh in atom.GetNeighbors() if self._is_count_valid(neigh)]
            n_neighs = len(neighs)
            if atom.GetAtomicNum() > 8:  # next row.
                # weird chemistry... likely wrong!
                penalties[i,:] = 2.
                penalties[:,i] = 2.
            elif n_neighs == 2:
                pass # no penalty!
            elif atom.GetIsAromatic():
                penalties[i, :] = 2 # this would result in a ring downgrade...
                penalties[:, i] = 2
            elif n_neighs == 3:
                penalties[i, :] = 1.  # 4 bonded carbon is not nice...
                penalties[:, i] = 1.
            else:
                penalties[i, :] = np.nan  # this will likely crash things.
                penalties[:, i] = np.nan
        # get closest pair.
        pendist_matrix = penalties + distance_matrix
        pendistance = np.nanmin(pendist_matrix)
        if np.isnan(pendistance):
            raise ValueError('This is impossible...')
        else:  # bonded
            p = np.where(pendist_matrix == pendistance)
            a = int(p[0][0])
            b = int(p[1][0])
            # absorb or bond
            distance = distance_matrix[a, b]
            penalty = penalties[a, b]
            if distance > 4:
                raise ValueError(f'The bond between {a} and {b} too long {distance} '+\
                                 f'from {indices_ring} and {indices_other}')
            elif distance > absorption_distance:
                # get bond type
                present_bond = mol.GetBondBetweenAtoms(ringcore.GetIdx(), other.GetIdx())
                bt = present_bond.GetBondType()
                if bt is None or bt == Chem.BondType.UNSPECIFIED:
                    bt = Chem.BondType.SINGLE
                mol.AddBond(a, b, bt)
                new_bond = mol.GetBondBetweenAtoms(a, b)
                BondProvenance.set_bond(new_bond, 'original')
                BondProvenance.copy_bond(present_bond, new_bond)
                # This is no longer required:
                # atom_a = mol.GetAtomWithIdx(a)
                # atom_b = mol.GetAtomWithIdx(b)
                #self._add_bond_if_possible(mol, atom_a, atom_b)
                self.journal.info(f'A novel bonding to ring was added {distance} {penalty}')
            else:
                # absorb the non-ring atom!
                self.journal.info(f'An atom was absorbed to ring was added {distance} {penalty}')
                if mol.GetAtomWithIdx(a).GetIntProp('_ori_i') == -1:
                    self._copy_bonding(mol, a, b)
                    self._mark_for_deletion(mol, b)
                    return [(a, b)]
                else:
                    self._copy_bonding(mol, b, a)
                    self._mark_for_deletion(mol, a)
                    return [(b, a)]

    # ======== Emergency ===============================================================================================

    def _emergency_joining(self, mol):
        """
        The last check to see if the mol is connected, before being rectified (valence fixes).
        """
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        n = len(frags)
        if n == 1:
            return mol
        else:
            while n > 1:
                self.journal.warning(f'Molecule disconnected in {n} parts. Please inspect final product!')
                name = mol.GetProp('_Name')
                for i, frag in enumerate(frags):
                    frag.SetProp('_Name', f'name.{i}')
                # find which fragments are closest.
                #TODO use the distance_matrix = self._get_distance_matrix(..) code
                closeness = np.ones([n, n])
                closeness.fill(float('nan'))
                for a, b in itertools.combinations(list(range(n)), 2):
                    closeness[a, b] = self._find_closest(frags[a], frags[b])[3]
                p = np.where(closeness == np.nanmin(closeness))
                frags = list(frags)
                first = frags[p[0][0]]
                second = frags[p[1][0]]
                mol = self.join_neighboring_mols(first, second)
                frags.remove(first)
                frags.remove(second)
                for part in frags:
                    mol = Chem.CombineMols(mol, part)
                    mol.SetProp('_Name', name)
                frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
                n = len(frags)
            return Chem.RWMol(mol)


        #
        #
        #     distance = distance_matrix[anchor_A, anchor_B]
        #
        #
        # rpos = np.array(conf.GetAtomPosition(ringcore_A.GetIdx()))
        # npos = np.array(conf.GetAtomPosition(ringcore_B.GetIdx()))
        # if np.linalg.norm(rpos - npos) <= 4:
        #     pairs = []
        #     for G_idxs, ref_i in [(A_idxs, neigh.GetIdx()), (B_idxs, ringcore.GetIdx())]:
        #         tm = np.take(dm, G_idxs, 0)
        #         tm2 = np.take(tm, [ref_i], 1)
        #         p = np.where(tm2 == np.nanmin(tm2))
        #         f = G_idxs[int(p[0][0])]
        #         tm2[p[0][0], :] = np.ones(tm2.shape[1]) * float('inf')
        #         p = np.where(tm2 == np.nanmin(tm2))
        #         s = G_idxs[int(p[1][0])]
        #         pairs.append((f, s))
        #     # now determine which are closer
        #     if dm[pairs[0][0], pairs[1][0]] < dm[pairs[0][0], pairs[1][1]]:
        #         mergituri.append((pairs[0][0], pairs[1][0]))
        #         mergituri.append((pairs[0][1], pairs[1][1]))
        #     else:
        #         mergituri.append((pairs[0][0], pairs[1][1]))
        #         mergituri.append((pairs[0][1], pairs[1][0]))

    def _detriangulate(self, mol: Chem.RWMol) -> None:
        """
        Prevents novel cyclopropanes and cyclobutanes.

        :param mol:
        :return:
        """
        for atom in mol.GetAtoms():
            atom_i = atom.GetIdx()
            for neigh in atom.GetNeighbors():
                neigh_i = neigh.GetIdx()
                if neigh_i < atom_i: # dont check twice...
                    continue
                # de triangulate
                third_i = self._get_triangle(atom, neigh)
                if third_i is not None:
                    # it is a triangle
                    third = mol.GetAtomWithIdx(third_i)
                    self.journal.debug(f'Triangle present {(atom_i, neigh_i, third_i)}.')
                    self._detriangulate_inner(mol,
                                              atoms=[atom, neigh, third],
                                              atom_indices=[atom_i, neigh_i, third_i],
                                              combinator=partial(itertools.combinations, r=2)
                                              )
                # de square-ify
                sq = self._get_square(atom, neigh)
                if sq is not None:
                    far_i, close_i = sq
                    far = mol.GetAtomWithIdx(far_i)
                    close = mol.GetAtomWithIdx(close_i) # second neighbour of atom
                    # bonding is:
                    # atom - neigh - far - close - atom
                    self.journal.debug(f'Square present {(atom_i, neigh_i, far_i, close_i)}.')
                    # combinations would fail at a atom - far bond
                    # the order is irrelevant if the same fun is called
                    combinator = lambda l: [[l[0], l[1]], [l[0], l[2]], [l[1], l[3]], [l[2], l[3]]]
                    self._detriangulate_inner(mol,
                                              atoms=[atom, neigh, close, far],
                                              atom_indices=[atom_i, neigh_i, close_i, far_i],
                                              combinator=combinator
                                              )

    def _detriangulate_inner(self,
                             mol: Chem.RWMol,
                             atoms: List[Chem.Atom],
                             atom_indices: List[int],
                             combinator: Callable):
        """
        Triangle/square agnostic.

        :param mol:
        :param atoms:
        :param atom_indices:
        :param combinator:
        :return:
        """
        bonds = [mol.GetBondBetweenAtoms(a, b) for a, b in combinator(atom_indices)]
        if any([bond is None for bond in bonds]):
            self.journal.critical(f'IMPOSSIBLE ERROR: detriangulate missing bond. {atom_indices}, {atoms}')
            return None
        provenances = BondProvenance.get_bonds(bonds)
        # original
        originality = [p == BondProvenance.ORIGINAL for p in provenances]
        if sum(originality) == 3:
            self.journal.warning('Triangle from original present. Kept, unless rectifiers is set to not tolerate')
        # length
        BL = partial(Chem.rdMolTransforms.GetBondLength, conf=mol.GetConformer())
        lengths = [BL(iAtomId=b.GetBeginAtomIdx(), jAtomId=b.GetEndAtomIdx()) for b in bonds]
        scores = np.array(lengths)
        scores[originality] = np.nan
        # atom properties. overbonded, ring etc.
        for fun, weight in self.closeness_weights:
            funscore = [fun(a) for a in atoms]
            scores += np.sum(list(combinator(funscore)), axis=1)
        d = np.nanmax(scores)
        if np.isnan(d):
            return None # impossible but okay.
        doomed_i = int(np.where(scores == d)[0])
        doomed_bond = bonds[doomed_i]
        a, b = doomed_bond.GetBeginAtomIdx(), doomed_bond.GetEndAtomIdx()
        self.journal.debug(f'Removing triangle/square forming bond between {a} and {b}')
        mol.RemoveBond(a, b)







