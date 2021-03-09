########################################################################################################################

__doc__ = \
    """
This is the ring collapsing code.

In Pymol, a ring diameter is 2.7894375324249268 Å.

    fab F
    print cmd.distance('/obj01///PHE`1/CE1','/obj01///PHE`1/CD2')
    """

########################################################################################################################

import itertools
import json
from collections import defaultdict
from functools import partial
from typing import Optional, Dict, List, Any, Tuple, Union, Callable

import numpy as np
from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D

from ._join_neighboring import _MonsterJoinNeigh
from .bond_provenance import BondProvenance


########################################################################################################################

class _MonsterRing( _MonsterJoinNeigh):

    def collapse_mols(self, mols: List[Chem.Mol]):
        mols = [self.collapse_ring(mol) for mol in mols]
        [self.offset(mol) for mol in mols]
        return mols

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
        self.keep_copy(mol, 'Rings expanded and original bonding restored.')
        self._add_novel_bonding(mol, rings)  # formerly `_ring_overlap_scenario` and `_infer_bonding_by_proximity`.
        self._delete_collapsed(mol)
        self._detriangulate(mol)
        try:
            mol = self._emergency_joining(mol)  # does not modify in place!
        except ConnectionError as error:
            if self.throw_on_discard:
                raise error
            else:
                self.journal.info('Disconnect ignored due to keep_all=False')
                mol = self.get_largest_fragment(mol)
        if mol is None:
            raise ValueError('(Impossible) Failed at some point...')
        elif isinstance(mol, Chem.RWMol):
            return mol.GetMol()
        else:
            return mol

    # =========== Offset ===============================================================================================



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
                pass  # ringcores have -1 ori_i
        # sort the ringcore
        for atom in self._get_collapsed_atoms(mol):
            old = json.loads(atom.GetProp('_ori_is'))
            new = [i + self._collapsed_ring_offset for i in old]
            atom.SetProp('_ori_is', json.dumps(new))
            old2new = {**old2new, **dict(zip(old, new))}
        # this has to be done afterwards in case of a bonded mol
        for atom in self._get_collapsed_atoms(mol):
            old_neighss = json.loads(atom.GetProp('_neighbors'))  # if i in old2new else i
            new_neighss = [[old2new[i] for i in old_neighs if i in old2new] for old_neighs in old_neighss]
            atom.SetProp('_neighbors', json.dumps(new_neighss))
        # determine if the new atoms have close neighbours.
        pass

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

    def _get_expansion_for_atom(self, ring: Dict[str, List[Any]], i: int) -> Dict[str, Any]:
        """
        ``_get_expansion_data`` returns from a mol the "expansion data for the rings"
        ``_get_expansion_for_atom`` given one of the list of the data from the latter (representing a ring core)
        and an index of which of the internal atoms that were collapsed return a dictionary of details
        of that atom.

        :param ring: see ``_get_expansion_data``
        :param i: the internal index. Say 'elements': ['C', 'C', 'C', 'O', 'C', 'C'].  i = 3 would will be Oxygen.
        :return:
        """
        try:
            return {k.replace('s', ''): ring[k][i] if isinstance(ring[k], list) else ring[k] for k in ring}
        except IndexError:
            troublesome = [k for k in ring if isinstance(ring[k], list) and len(ring[k]) <= i]
            if len(troublesome) == 0:
                raise IndexError(f'There is a major issue with ring data for index {i}: {ring}')
            elif troublesome[0] == 'current_is':
                self.journal.warning(f'One atom lacks a current index!' + \
                                     'This is a fallback that should not happen')
                mol = ring['atom'].GetOwningMol()
                ring['current_is'] = [self._get_new_index(mol, old_i, search_collapsed=False) for old_i in
                                      ring['ori_is']]
                return self._get_expansion_for_atom(ring, i)
            else:
                raise IndexError(f'The indices of the collapsed atom do not extend to {i} for {troublesome}')

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
            ringcore = ring['atom']
            indices = []  # will store current indices
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
                    natom.SetIntProp('_ring_i', ringcore.GetIdx())
                    indices.append(n)
            ringcore.SetIntProp('_ring_i', ringcore.GetIdx())  # it really really should have not changed.
            ringcore.SetProp('_current_is', json.dumps(indices))
            ring['current_is'] = indices

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
        self.journal.debug('Adding novel bonding (if any)...')
        # ===== Deal with Ring on ring bonding ------------------------------------------
        novel_ringcore_pairs = self._get_novel_ringcore_pairs(mol, rings, cutoff=1.5)
        # these is a list of Chem.Atom pairs.
        for ringcore_A, ringcore_B in novel_ringcore_pairs:
            self.journal.debug('determining novel bond between ring markers ' + \
                               f'{ringcore_A.GetIdx()} and {ringcore_B.GetIdx()}')
            # _determine_mergers_novel_ringcore_pair finds mergers
            self._determine_mergers_novel_ringcore_pair(mol, ringcore_A, ringcore_B)
        # ===== Deal with Ring on other bonding ------------------------------------------
        # formerly: _infer_bonding_by_proximity
        novel_other_pairs = self._get_novel_other_pairs(mol, rings, 1.0)
        for ringcore, other in novel_other_pairs:
            self.journal.debug(f'determining novel bond between ' + \
                               f'ring marker {ringcore.GetIdx()} and non-ring {other.GetIdx()}')
            # _determine_mergers_novel_ringcore_pair finds, bonds and marks for deletion.
            self._determine_mergers_novel_other_pair(mol, ringcore, other)
        # ===== Clean up ------------------------------------------------------------------
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

    def _get_novel_ringcore_pairs(self,
                                  mol: Chem.Mol,
                                  rings: List[Dict[str, List[Any]]],
                                  cutoff: float) \
            -> List[Tuple[Chem.Atom, Chem.Atom]]:
        """
        Get the ring atoms that are bonded or closer than a cutoff.
        It is the countepart to _get_novel_ringcore_pairs

        :param mol:
        :param rings: output of `_get_expansion_data`. See that for details.
        :type rings: List[Dict[str, List[Any]]]
        :param cutoff:
        :return:
        """
        # ----------------------------------------------------------
        # scenario where they are closer than the cutoff
        close_pairs = self._get_close_novel_ringcores(mol, rings, cutoff)  # 2.7889 ring diameter + 1.45 C-C bond
        # ----------------------------------------------------------
        # scenario where they are bonded...
        bonded_pairs = self._get_novel_ringcore_bonded_pairs(rings)
        # ------------ merge lists -----------
        return self.merge_pairing_lists(close_pairs + bonded_pairs)

    def _get_novel_ringcore_bonded_pairs(self, rings):
        """
        called by _get_novel_ringcore_pairs alongside _get_close_novel_ringcores
        Opposite of _get_novel_other_bonded_pairs

        :param rings:
        :return:
        """
        bonded_pairs = []
        for ring in rings:
            # ring: Dict[str, List[Any]]
            ringcore = ring['atom']  # Chem.Atom the ring core marker, not a real atom.
            for neigh in ringcore.GetNeighbors():
                # find those ringcore atoms that are connected. via these conditions:
                has_ringcore_neighbor = neigh.HasProp('_ori_name') and neigh.GetIntProp('_ori_i') == -1  # -1 is ring.
                is_new_pair = ringcore.GetIdx() < neigh.GetIdx()  # to avoid doing it twice each direction
                is_novel_connection = neigh.HasProp('_ori_name') and ring['ori_name'] != neigh.GetProp('_ori_name')
                # checking:
                if has_ringcore_neighbor and is_new_pair and is_novel_connection:
                    # This ringcore atom shares a novel border with another ringcore atom
                    bonded_pairs.append((ringcore, neigh))  # i.e. List[Tuple[Chem.Atom, Chem.Atom]}
        return bonded_pairs

    def merge_pairing_lists(self,
                            nonunique_pairs: List[Tuple[Chem.Atom, Chem.Atom]],
                            ringcore_first=True) \
            -> List[Tuple[Chem.Atom, Chem.Atom]]:
        # complicate because mol.GetAtomWithIdx(2) == mol.GetAtomWithIdx(2) is False
        seen = []
        pairs = []
        for atom_a, atom_b in nonunique_pairs:
            # sort out
            ai = atom_a.GetIdx()
            bi = atom_b.GetIdx()
            if ai == bi:
                self.journal.debug(f'Bond to self incident with {ai} (ring? {atom_a.HasProp("_current_is") == 1})')
                continue
            if ringcore_first and atom_a.HasProp('_current_is') and not atom_b.HasProp('_current_is'):
                ringcore = ai
                other = bi
                pairing = f'{ringcore}-{other}'
            elif ringcore_first and atom_b.HasProp('_current_is') and not atom_a.HasProp('_current_is'):
                ringcore = bi
                other = ai
                pairing = f'{ringcore}-{other}'
                atom_a, atom_b = atom_b, atom_a
            elif atom_a.HasProp('_current_is') and atom_b.HasProp('_current_is'):
                low_i, high_i = sorted([ai, bi])
                pairing = f'{low_i}-{high_i}'
            else:
                # these should not have been let thorugh but other novel does not filter them.
                # self.journal.debug(f'Non-ring to non-ring closeness flagged!')
                continue
            # verify unseen
            if pairing in seen:
                pass
            else:
                seen.append(pairing)
                pairs.append((atom_a, atom_b))
        return pairs

    def _get_ring_atom_indices_per_origin(self, rings: List[Dict[str, List[Any]]]) -> Dict[str, List[int]]:
        atomdex = defaultdict(set)
        for ring in rings:
            origin_name = ring['ori_name']  # ring['ori_name'] is same as ringcore.GetProp('_ori_name')
            atomdex[origin_name].update(ring['current_is'])  # ring['current_is'] = json ringcore.GetProp('_current_is')
        return atomdex

    def _get_atom_indices_per_origin(self, mol: Chem.Mol) -> Dict[str, List[int]]:
        atomdex = defaultdict(list)
        for atom in mol.GetAtoms():
            if not atom.HasProp('ori_name'):
                atomdex['unknown'].append(atom.GetIdx())
            else:
                name = atom.GetProp('ori_name')
                atomdex[name].append(atom.GetIdx())
        return atomdex

    def _get_close_novel_ringcores(self, mol: Chem.Mol, rings: List[Dict[str, List[Any]]], cutoff: float):
        cnrai = self._get_close_novel_ring_atoms_indices(mol, rings, cutoff)
        return self._indices_to_atoms_n_cores(mol, cnrai)

    def _indices_to_atoms_n_cores(self, mol: Chem.Mol, index_pairs: List[Tuple[int, int]]):
        """
        Give a list of index pairs convert them to a list of pairs of atoms/cores
        """
        # these are indices of ring atoms, not ring cores
        get_atom = lambda i: mol.GetAtomWithIdx(int(i))
        get_ringcore = lambda atom: get_atom(atom.GetIntProp('_ring_i')) if atom.HasProp('_ring_i') else atom
        idx2ringcore = lambda i: get_ringcore(get_atom(i))
        return [(idx2ringcore(ai), idx2ringcore(bi)) for ai, bi in index_pairs]

    def _get_close_novel_ring_atoms_indices(self,
                                            mol: Chem.Mol,
                                            rings: List[Dict[str, List[Any]]],
                                            cutoff: int) -> List[Tuple[int, int]]:
        """
        Get the list of pairs of indices derived from a ring that are closer that ``cutoff`` to another ring atom.
        Note, the operations are between real ring atoms not ring core markers.
        Hence why ``_get_close_novel_ringcores`` does a conversion.


        :param mol:
        :param rings: output of `_get_expansion_data`. See that for details.
        :param cutoff:
        :return:
        """
        atomdex = self._get_ring_atom_indices_per_origin(rings)
        # not calling get_distance matrxi because tehre may be more than 2 origins.
        distance_matrix = Chem.Get3DDistanceMatrix(mol)
        for origin_name, indices in atomdex.items():
            # this blanks subsquares of same origin but not interesections of different indices from atomdex...
            self._nan_fill_submatrix(distance_matrix, list(indices))
        self._nan_fill_others(mol, distance_matrix, [idx for idcs in atomdex.values() for idx in idcs])
        return self._get_closest_from_matrix(distance_matrix, cutoff)

    def _get_closest_from_matrix(self, matrix: np.ndarray, cutoff: int) -> List[Tuple[int, int]]:
        # get the pair of atom indices that are less thna cutoff.
        # where returns a tuple of np.arrays of dtype=np.int64
        with np.errstate(invalid='ignore'):
            return list(zip(*[w.astype(int) for w in np.where(matrix < cutoff)]))

    def _get_close_novel_ring_other_indices(self,
                                            mol: Chem.Mol,
                                            rings: List[Dict[str, List[Any]]],
                                            cutoff: int) -> List[Tuple[int, int]]:
        """
        Get the list of pairs of indices between an ring atom and a non-ring atom from a different origin that is too close.

        :param mol:
        :param rings: output of `_get_expansion_data`. See that for details.
        :param cutoff:
        :return:
        """
        atomdex = self._get_ring_atom_indices_per_origin(rings)
        oridex = self._get_atom_indices_per_origin(mol)
        distance_matrix = Chem.Get3DDistanceMatrix(mol)
        # blank all rings to rings
        self._nan_fill_submatrix(distance_matrix, [i for l in atomdex.values() for i in l])
        closeness = []
        for origin_name, indices in atomdex.items():
            sub = distance_matrix.copy()
            self._nan_fill_submatrix(sub, oridex[origin_name])
            # get the pair of atom indices that are less thna cutoff.
            # where returns a tuple of np.arrays of dtype=np.int64
            closeness.extend(self._get_closest_from_matrix(distance_matrix, cutoff))
        return closeness

    def _determine_mergers_novel_ringcore_pair(self,
                                               mol: Chem.RWMol,
                                               ringcore_A: Chem.Atom,
                                               ringcore_B: Chem.Atom) -> List[Tuple[int, int]]:
        """
        Preps to resolve ringcore pairs.
        Formerly part of ``_ring_overlap_scenario``.
        Preps without deleting. ``_delete_marked(mol)`` does that.
        bonded, spiro, fused

        :param mol:
        :param ringcore_A:
        :param ringcore_B:
        :return: list of atoms to be merged
        """
        absorption_distance = 1.  # Å
        # print('A', ringcore_A, ringcore_A.GetIdx(), ringcore_A.GetIntProp('_ori_i'))
        # print('B', ringcore_B, ringcore_B.GetIdx(), ringcore_B.GetIntProp('_ori_i'))
        indices_A = json.loads(ringcore_A.GetProp('_current_is'))
        indices_B = json.loads(ringcore_B.GetProp('_current_is'))
        distance_matrix = self._get_distance_matrix(mol, indices_A, indices_B)  # currently in `_join_neighboring`.
        # distance matrix is for the whole thing
        # TODO merge into _get_distance_matrix
        # blanking the other atom indices:
        self._nan_fill_others(mol, distance_matrix, indices_A + indices_B)
        # get closest pair.
        distance = np.nanmin(distance_matrix)
        if np.isnan(distance):
            self.journal.critical('This is impossible. Two neighbouring rings cannot be connected.')
            return []
        elif distance > absorption_distance:  # bonded
            p = np.where(distance_matrix == distance)
            a = int(p[0][0])
            b = int(p[1][0])
            present_bond = mol.GetBondBetweenAtoms(ringcore_A.GetIdx(), ringcore_B.GetIdx())
            self._add_bond_by_reference(mol, a, b, present_bond)
            self.journal.info('A novel bond-connected ring pair was found')
            self._mark_for_deletion(mol, b)
            self._copy_bonding(mol, a, b, force=True)
            return []  # bonded
        else:  # Spiro or fused.
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
                return [(a, b), (c, d)]  # fused
            else:
                self.journal.info('A novel spiro ring pair was found')
                return [(a, b)]  # spiro

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

    def _get_novel_other_bonded_pairs(self, rings) -> List[Tuple[Chem.Atom, Chem.Atom]]:
        pairs = []
        for ring in rings:
            # ring: Dict[str, List[Any]]
            ringcore = ring['atom']  # Chem.Atom
            # --------- find ring atom - non-origin other pairs that are bonded
            for neigh in ringcore.GetNeighbors():
                # find those ringcore-other atoms that are connected. via these conditions:
                has_notringcore_neighbor = (not neigh.HasProp('_ori_name')) or neigh.GetIntProp('_ori_i') != -1
                is_novel_connection = (not neigh.HasProp('_ori_name')) or ring['ori_name'] != neigh.GetProp('_ori_name')
                if has_notringcore_neighbor and is_novel_connection:
                    # This ringcore atom shares a novel border with another ringcore atom
                    pairs.append((ringcore, neigh))
        return pairs

    def _get_novel_other_pairs(self, mol, rings, cutoff: float) -> List[Tuple[Chem.Atom, Chem.Atom]]:
        """
        similar to self._get_novel_ringcore_pairs... but opposite
        It deals with bonded and close pairs.

        :param mol:
        :param rings: output of `_get_expansion_data`. See that for details.
        :type rings: List[Dict[str, List[Any]]]
        :param cutoff:
        :return:
        """
        # ----------------------------------------------------------
        # scenario where they are closer than the cutoff
        close_pairs = self._get_close_novel_others(mol, rings, cutoff)  # 2.7889 ring diameter + 1.45 C-C bond
        # ----------------------------------------------------------
        # scenario where they are bonded...
        bonded_pairs = self._get_novel_other_bonded_pairs(rings)
        # ------------ merge lists -----------
        return self.merge_pairing_lists(close_pairs + bonded_pairs)

    def _get_close_novel_others(self,
                                mol,
                                rings,
                                cutoff):
        idx_pairs = self._get_close_novel_ring_other_indices(mol, rings, cutoff)
        # filter out those that are bonded already to core atom.
        return self._indices_to_atoms_n_cores(mol, idx_pairs)

    def _get_merging_penalties(self, mol, shape: Tuple[int, int], indices_ring):
        """
        confusingly this is a different set of penalties to `_get_joining_penalties` which is for joining
        :param mol:
        :param shape:
        :param indices_ring: The indices inside the ring atom, aka. prop _current_is
        :return:
        """
        penalties = np.zeros(shape)
        for i in indices_ring:
            atom = mol.GetAtomWithIdx(i)
            neighs = [neigh for neigh in atom.GetNeighbors() if self._is_count_valid(neigh)]
            n_neighs = len(neighs)
            if atom.GetAtomicNum() > 8:  # next row.
                # weird chemistry... likely wrong!
                penalties[i, :] = 2.
                penalties[:, i] = 2.
            elif n_neighs == 2:
                pass  # no penalty!
            elif atom.GetIsAromatic():
                penalties[i, :] = 2  # this would result in a ring downgrade...
                penalties[:, i] = 2
            elif n_neighs == 3:
                penalties[i, :] = 1.  # 4 bonded carbon is not nice...
                penalties[:, i] = 1.
            else:
                penalties[i, :] = np.nan  # this will likely crash things.
                penalties[:, i] = np.nan
        return penalties

    def _get_distance(self, atom_a: Chem.Atom, atom_b: Chem.Atom) -> np.float:
        """
        Not sure where doing it manually is quicker than getting the whole 3D distance table.

        :param atom_a:
        :param atom_b:
        :return:
        """
        conf = atom_a.GetOwningMol().GetConformer()
        get_pos = lambda atom: np.array(conf.GetAtomPosition(atom.GetIdx()))
        return np.linalg.norm(get_pos(atom_a) - get_pos(atom_b))


    def _determine_mergers_novel_other_pair(self,
                                            mol: Chem.RWMol,
                                            ringcore: Chem.Atom,
                                            other: Chem.Atom) -> List[Tuple[int, int]]:
        """
        Like _determine_mergers_novel_ringcore_pair finds, bonds and marks for deletion.
        It however finds atoms to absorb between a ring and a given non-ring atom.

        :param mol:
        :param ringcore:
        :param other:
        :return:
        """
        # ---- Prep data.
        indices_ring = json.loads(ringcore.GetProp('_current_is'))
        index_other = other.GetIdx()
        indices_other = [index_other]
        index_core = ringcore.GetIdx()
        distance_matrix = self._get_distance_matrix(mol, indices_ring,
                                                    indices_other)  # currently in `_join_neighboring`.
        self._nan_fill_others(mol, distance_matrix, indices_ring + indices_other)
        # merging penalties
        penalties = self._get_merging_penalties(mol, distance_matrix.shape, indices_ring)
        core_absorption_distance = 1.5  # Å between ring core and other. 2.8 Å is diameter.
        core_other_distance = self._get_distance(ringcore, other)
        if core_other_distance < core_absorption_distance:
            # ------ Within ring must go. No penalties.
            self.journal.debug(f'(DetMergeNovOther *{index_core}, {index_other}). Other is within ring. Forcing absoption')
            absorption_distance = 9999
            pendist_matrix = distance_matrix
        else:
            # ------ Assess cases normally.
            absorption_distance = 1.  # Å between ring atom and other.
            pendist_matrix = penalties + distance_matrix
        # get closest pair.
        pendistance = np.nanmin(pendist_matrix)
        if np.isnan(pendistance):
            self.journal.warning(f'(DetMergeNovOther*{index_core}, {index_other}). This is impossible...')
            return []
        else:  # bonded
            p = np.where(pendist_matrix == pendistance)
            a = int(p[0][0])
            b = int(p[1][0])
            assert index_other in (a, b), 'CRITICIAL: Matrix error!'
            # absorb or bond
            distance = distance_matrix[a, b]
            penalty = penalties[a, b]  # penalties were already applied. this is for msgs only
            if distance > 4:
                self.journal.warning(f'(DetMergeNovOther: {a}, {b}). ' + \
                                     f'The bond between {a} and {b} too long {distance} ' + \
                                     f'from {indices_ring} and {indices_other}')
                return []
            elif distance > absorption_distance:
                self.journal.info(f'(DetMergeNovOther: {a}, {b}). A novel bonding to ring may be added. ' + \
                                  f'd: {distance} p: {penalty}')
                # get bond type
                present_bond = mol.GetBondBetweenAtoms(ringcore.GetIdx(), other.GetIdx())
                self._add_bond_by_reference(mol, a, b, present_bond)
            else:
                # absorb the non-ring atom!
                self.journal.info(f'(DetMergeNovOther: {a}, {b}). An atom was absorbed to ring. ' + \
                                  f'd: {distance} p: {penalty}')
                if mol.GetAtomWithIdx(a).GetIntProp('_ori_i') == -1:
                    self._copy_bonding(mol, a, b)
                    self._mark_for_deletion(mol, b)
                    return [(a, b)]
                else:
                    self._copy_bonding(mol, b, a)
                    self._mark_for_deletion(mol, a)
                    return [(b, a)]

    def _add_bond_by_reference(self, mol, a, b, reference_bond):
        """
        _copy_bonding copies all the bonds. THis just adds one like the reference bond.
        It calls ``_add_bond_if_possible`` if its a closeness bond, i.e. reference_bond is None
        It calls ``_add_bond_regardlessly`` if its an orginal one.

        :param mol:
        :param a:
        :param b:
        :param reference_bond:
        :return:
        """
        atom_a = mol.GetAtomWithIdx(a)
        atom_b = mol.GetAtomWithIdx(b)
        if reference_bond is None:
            self._add_bond_if_possible(mol, atom_a, atom_b, 'other_novel')
        else:
            bt = reference_bond.GetBondType()
            if bt is None or bt == Chem.BondType.UNSPECIFIED:
                bt = Chem.BondType.SINGLE
            self._add_bond_regardlessly(mol, atom_a, atom_b, bt, BondProvenance.get_bond(reference_bond).name)

    # ======== Emergency ===============================================================================================

    def _emergency_joining(self, mol: Chem.Mol) -> Chem.Mol:
        return self._join_internally(mol, severe=True)

    def _join_internally(self, mol: Chem.Mol, severe: bool=False) -> Chem.Mol:
        """
        The last check to see if the mol is connected.
        This differs (and calls) ``join_neighboring_mols``

        """
        is_rw = isinstance(mol, Chem.RWMol)
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        n = len(frags)
        if n == 1:
            return mol
        else:
            while n > 1:
                if severe:
                    self.journal.warning(f'Molecule disconnected in {n} parts. Please inspect final product!')
                else:
                    self.journal.debug('Linking two disconnected fragments')
                # ----- get names ---------------
                name = mol.GetProp('_Name')
                for i, frag in enumerate(frags):
                    frag.SetProp('_Name', f'name.{i}')
                # find which fragments are closest ------------------------------
                # TODO use the distance_matrix = self._get_distance_matrix(..) code
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
                # ---- reset variables ---------
                for part in frags:
                    mol = Chem.CombineMols(mol, part)
                    mol.SetProp('_Name', name)
                frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
                n = len(frags)
            if is_rw:
                return Chem.RWMol(mol)
            else:
                return mol

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
                if neigh_i < atom_i:  # dont check twice...
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
                    close = mol.GetAtomWithIdx(close_i)  # second neighbour of atom
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
            return None  # impossible but okay.
        doomed_i = int(np.where(scores == d)[0])
        doomed_bond = bonds[doomed_i]
        a, b = doomed_bond.GetBeginAtomIdx(), doomed_bond.GetEndAtomIdx()
        self.journal.debug(f'Removing triangle/square forming bond between {a} and {b}')
        mol.RemoveBond(a, b)
