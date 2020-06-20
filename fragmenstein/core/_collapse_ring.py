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
from typing import Optional, Dict, List, Any, Tuple
import numpy as np
from collections import Counter

import logging

log = logging.getLogger('Fragmenstein')


class Ring:
    def __init__(self, _debug_draw=False):
        self._debug_draw = _debug_draw

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

    # =========== Offset ===============================================================================================

    _collapsed_ring_offset = 0

    def _offset_collapsed_ring(self, mol: Chem.Mol):
        """
        This is to prevent clashes.
        The numbers of the ori indices stored in collapsed rings are offset by the class variable (_collapsed_ring_offset)
        multiples of 100. (autoincrements to avoid dramas)

        :param mol:
        :return:
        """
        self._collapsed_ring_offset += 100
        for atom in self._get_collapsed_atoms(mol):
            old = json.loads(atom.GetProp('_ori_is'))
            new = [i + self._collapsed_ring_offset for i in old]
            atom.SetProp('_ori_is', json.dumps(new))
            old2new = dict(zip(old, new))
            if self._debug_draw:
                print('UPDATE', old2new)
            old_neighss = json.loads(atom.GetProp('_neighbors'))
            new_neighss = [[old2new[i] if i in old2new else i for i in old_neighs] for old_neighs in old_neighss]
            atom.SetProp('_neighbors', json.dumps(new_neighss))

    def _offset_origins(self, mol: Chem.Mol):
        """
        This is to prevent clashes.

        :param mol:
        :return:
        """
        self._collapsed_ring_offset += 100
        old2new = {}
        for atom in mol.GetAtoms():
            if atom.GetIntProp('_ori_i') != -1:
                o = atom.GetIntProp('_ori_i')
                n = o + self._collapsed_ring_offset
                atom.SetIntProp('_ori_i', n)
                old2new[o] = n
        for atom in self._get_collapsed_atoms(mol):
            old = json.loads(atom.GetProp('_ori_is'))
            new = [i + self._collapsed_ring_offset for i in old]
            atom.SetProp('_ori_is', json.dumps(new))
            old2new = {**old2new, **dict(zip(old, new))}
            if self._debug_draw:
                print('UPDATE', old2new)
            old_neighss = json.loads(atom.GetProp('_neighbors'))
            new_neighss = [[old2new[i] if i in old2new else i for i in old_neighs] for old_neighs in old_neighss]
            atom.SetProp('_neighbors', json.dumps(new_neighss))

    def _renumber_original_indices(self, mol: Chem.Mol,
                                   mapping: Dict[int, int],
                                   name_restriction: Optional[str] = None):
        """

        :param mol:
        :param mapping: old index to new index dictionary
        :param name_restriction:
        :return:
        """
        for atom in mol.GetAtoms():
            if name_restriction is not None and atom.HasProp('_ori_name') and atom.GetProp(
                    '_ori_name') != name_restriction:
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

    # =========== Exand ================================================================================================

    def _get_expansion_data(self, mol: Chem.Mol) -> List[Dict[str, List[Any]]]:
        """
        Returns a list for each collapsed ring marking atom each with a dictionary

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

    def expand_ring(self, mol: Chem.Mol, bonded_as_original=False):
        """
        Undoes collapse ring

        :param mol:
        :param bonded_as_original: if false bonds by proximity, if true bonds as remembered
        :return:
        """
        mol = Chem.RWMol(mol)
        rings = self._get_expansion_data(mol)
        self._place_ring_atoms(mol, rings)
        mergituri = self._ring_overlap_scenario(mol, rings)
        # bonding
        if bonded_as_original:
            self._restore_original_bonding(mol, rings)
            self._fix_overlap(mol, mergituri)
            self._delete_collapsed(mol)
        else:
            self._connenct_ring(mol, rings)
            self._mark_neighbors(mol, rings)
            self._fix_overlap(mol, mergituri)
            self._delete_collapsed(mol)
            self._infer_bonding_by_proximity(mol) # absorb_overclose and join_overclose
        # this should not happen... but it can!
        mol = self._emergency_joining(mol)
        # _emergency_joining returns a Chem.Mol not a Chem.RWMol
        # prevent weird nested rings.
        mol = self._prevent_conjoined_ring(mol)
        mol = self._prevent_weird_rings(mol)
        if mol is None:
            raise ValueError('(Impossible) Failed at some point...')
        elif isinstance(mol, Chem.RWMol):
            return mol.GetMol()
        else:
            return mol

    def _emergency_joining(self, mol):
        """
        The last check to see if the mol is connected, before being rectified (valence fixes).
        """
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        n = len(frags)
        while n > 1:
            log.warning(f'Molecule disconnected in {n} parts. Please inspect.')
            name = mol.GetProp('_Name')
            for i, frag in enumerate(frags):
                frag.SetProp('_Name', f'name.{i}')
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
        return mol

    def _fix_overlap(self, mol, mergituri):
        morituri = []
        for pair in mergituri:
            self._absorb(mol, *pair)
            morituri.append(pair[1])
        for i in sorted(morituri, reverse=True):
            mol.RemoveAtom(i)

    def _prevent_conjoined_ring(self, mol: Chem.Mol) -> Chem.Mol:
        """
        This kills bridging bonds with not atoms in the bridge within rings.
        So it is bridged, fused and spiro safe.
        It removes only one bond, so andamantane/norbornane are safe.
        :param mol:
        :return:
        """
        c = Counter([i for ring in self._get_ring_info(mol) for i in ring])
        nested = [k for k in c if c[k] >= 3]
        pairs = [(idx_a, idx_b) for idx_a, idx_b in itertools.combinations(nested, r=2) if
                 mol.GetBondBetweenAtoms(idx_a, idx_b) is not None]
        rank = sorted(pairs, key=lambda x: c[x[0]] + c[x[1]], reverse=True)
        if len(rank) > 0:
            idx_a, idx_b = rank[0]
            if not isinstance(mol, Chem.RWMol):
                mol = Chem.RWMol(mol)
            mol.RemoveBond(idx_a, idx_b)
            log.info(f'Zero-atom bridged ring issue: bond between {idx_a}-{idx_b} removed')
            return self._prevent_conjoined_ring(mol)
        elif isinstance(mol, Chem.RWMol):
            return mol.GetMol()
        else:
            return mol

    def _place_between(self, mol: Chem.RWMol, a: int, b: int, aromatic=True):
        oribond = mol.GetBondBetweenAtoms(a, b)
        if oribond is None:
            print('FAIL')
            return None  # fail
        elif aromatic:
            bt = Chem.BondType.AROMATIC
        else:
            bt = oribond.GetBondType()
        idx = mol.AddAtom(Chem.Atom(6))
        neoatom = mol.GetAtomWithIdx(idx)
        atom_a = mol.GetAtomWithIdx(a)
        atom_b = mol.GetAtomWithIdx(b)
        if aromatic:
            neoatom.SetIsAromatic(True)
            atom_a.SetIsAromatic(True)
            atom_b.SetIsAromatic(True)
        # prevent constraints
        neoatom.SetBoolProp('_Novel', True)
        atom_a.SetBoolProp('_Novel', True)
        atom_b.SetBoolProp('_Novel', True)
        # fix position
        conf = mol.GetConformer()
        pos_A = conf.GetAtomPosition(a)
        pos_B = conf.GetAtomPosition(b)
        x = pos_A.x / 2 + pos_B.x / 2
        y = pos_A.y / 2 + pos_B.y / 2
        z = pos_A.z / 2 + pos_B.z / 2
        conf.SetAtomPosition(idx, Point3D(x, y, z))
        # fix bonds
        mol.RemoveBond(a, b)
        mol.AddBond(a, idx, bt)
        mol.AddBond(b, idx, bt)

    def _prevent_weird_rings(self, mol: Chem.Mol):
        if not isinstance(mol, Chem.RWMol):
            mol = Chem.RWMol(mol)
        ringatoms = self._get_ring_info(mol) #GetRingInfo().AtomRings()
        for ring_A, ring_B in itertools.combinations(ringatoms, r=2):
            shared = set(ring_A).intersection(set(ring_B))
            if len(shared) == 0:
                log.debug('This molecule separate rings')
                pass  # separate rings
            elif len(shared) == 1:
                log.debug('This molecule has a spiro bicycle')
                pass  # spiro ring.
            elif len(shared) == 2:
                log.debug('This molecule has a fused ring')
                if mol.GetBondBetweenAtoms(*shared) is not None:
                    pass  # indole/naphtalene
                    small, big = sorted([ring_A, ring_B], key=lambda ring: len(ring))
                    if len(small) == 4:
                        log.warning('This molecule has a benzo-azetine–kind-of-thing: expanding to indole')
                        # Chem.MolFromSmiles('C12CCCCC1CC2')
                        # benzo-azetine is likely an error: add and extra atom
                        a, b = set(small).difference(big)
                        self._place_between(mol, a, b)
                    elif len(small) == 3:
                        log.warning('This molecule has a benzo-cyclopropane–kind-of-thing: expanding to indole')
                        # Chem.MolFromSmiles('C12CCCCC1C2')
                        # benzo-cyclopronane is actually impossible at this stage.
                        a = list(set(small).difference(big))[0]
                        for b in shared:
                            self._place_between(mol, a, b)
                    else:
                        pass  # indole and nathalene
                elif (len(ring_A), len(ring_B)) == (6, 6):
                    raise Exception('This is utterly impossible')
                else:
                    print(f'mysterious ring system {len(ring_A)} + {len(ring_B)}')
                    pass  # ????
            elif len(shared) < self.atoms_in_bridge_cutoff:
                #adamantene/norbornane/tropinone kind of thing
                log.warning('This molecule has a bridge: leaving')
                pass  # ideally check if planar...
            else:
                log.warning('This molecule has a bridge that will be removed')
                mol = self._prevent_bridge_ring(mol, ring_A)
                # start from scratch.
                return self._prevent_weird_rings(mol)
        return mol.GetMol()

    def _prevent_bridge_ring(self, mol: Chem.RWMol, examplar: Tuple[int]):
        ## This is really
        # examplar is ring
        ringatoms = self._get_ring_info(mol) #GetRingInfo().AtomRings()
        ringatoms = [ring for ring in ringatoms if set(ring).intersection(examplar)]
        ring_idx = list(range(len(ringatoms)))
        shared_count = {}
        for ra, rb in itertools.combinations(ring_idx, r=2):
            shared_count[(ra, rb)] = len(set(ringatoms[ra]).intersection(set(ringatoms[rb])))
        if len(shared_count) == 0:
            return mol
        ra, rb = list(shared_count.keys())[0]
        shared = list(set(ringatoms[ra]).intersection(ringatoms[rb]))
        pairs = [(a, b) for a, b in itertools.combinations(shared, r=2) if mol.GetBondBetweenAtoms(a, b) is not None]
        c = Counter([i for pair in pairs for i in pair])
        ring_A, ring_B = ringatoms[ra], ringatoms[rb]
        small, big = sorted([ring_A, ring_B], key=lambda ring: len(ring))
        inners = [i for i in c if c[i] > 1]
        x = list(set(shared).difference(inners))
        if x != 2:
            log.CRITICAL(f'This is impossible. {ringatoms} share {shared} with {inners} in the inside and {x} on the edge?')
            return mol
        a, b = x
        if len(big) > 6:
            log.warning(f'Removing {len(inners)} bridging atoms and replacing with fused ring')
            # bond the vertices
            bt = Chem.BondType.SINGLE # ???
            if mol.GetBondBetweenAtoms(a, b) is None:
                mol.AddBond(a, b, bt)
            else:
                log.warning('This is really odd! Why is there a bond already??')
            # remove the middle atoms.
            for i in sorted(inners, reverse=True):
                mol.RemoveAtom(i)
        else:
            log.warning(f'Shriking the smaller ring to change from bridged to fused.')
            # get the neighbour in the small atom to a vertex.
            neighs = [neigh for neigh in mol.GetAtomWithIdx(a).GetNeighbors() if
                      neigh.GetIdx() not in shared and neigh.GetIdx() in small]
            neigh = sorted(neighs, key=lambda atom: atom.GetSymbol() != 'C')[0]
            bt = mol.GetBondBetweenAtoms(a, neigh.GetIdx()).GetBondType()
            mol.RemoveBond(a, neigh.GetIdx())
            new_neigh = [neigh for neigh in mol.GetAtomWithIdx(a).GetNeighbors() if neigh.GetIdx() in shared][0]
            mol.AddBond(neigh.GetIdx(), new_neigh.GetIdx(), bt)
            neigh.SetBoolProp('_Novel', True)
            new_neigh.SetBoolProp('_Novel', True)
            mol.GetAtomWithIdx(a).SetBoolProp('_Novel', True)
        return mol


    def _place_ring_atoms(self, mol, rings):
        conf = mol.GetConformer()
        for ring in rings:
            # atom addition
            for i in range(len(ring['elements'])):
                d = self._get_expansion_for_atom(ring, i)
                if self._is_present(mol, d['ori_i']):
                    natom = self._get_new_index(mol, d['ori_i'], search_collapsed=False)
                    if self._debug_draw:
                        print(f"{natom} (formerly {d['ori_i']} existed already!")
                else:
                    n = mol.AddAtom(Chem.Atom(d['element']))
                    natom = mol.GetAtomWithIdx(n)
                    conf.SetAtomPosition(n, Point3D(d['x'], d['y'], d['z']))
                    natom.SetIntProp('_ori_i', d['ori_i'])
                    natom.SetDoubleProp('_x', d['x'])
                    natom.SetDoubleProp('_y', d['y'])
                    natom.SetDoubleProp('_z', d['z'])
                    natom.SetProp('_ori_name', d['ori_name'])

    def _ring_overlap_scenario(self, mol, rings):
        # resolve the case where a border of two rings is lost.
        # the atoms have to be ajecent.
        dm = Chem.Get3DDistanceMatrix(mol)
        mergituri = []
        for ring in rings:
            for n in ring['atom'].GetNeighbors():
                if n.GetIntProp('_ori_i') == -1 and ring['atom'].GetIdx() < n.GetIdx():  # it may share a border.
                    # variables to assess if overlap or bond
                    conf = mol.GetConformer()
                    rpos = np.array(conf.GetAtomPosition(ring['atom'].GetIdx()))
                    npos = np.array(conf.GetAtomPosition(n.GetIdx()))
                    if np.linalg.norm(rpos - npos) > 4: # is bond
                        # is it connected via a bond and not an overlapping atom?
                        # 2.8 ring dia vs. 2.8 + 1.5 CC bond
                        # this will be fixed depending on if from history or not.
                        pass
                    else: # is overlap
                        A_idxs_old = ring['ori_is']
                        A_idxs = [self._get_new_index(mol, i, search_collapsed=False) for i in A_idxs_old]
                        B_idxs_old = json.loads(n.GetProp('_ori_is'))
                        B_idxs = [self._get_new_index(mol, i, search_collapsed=False) for i in B_idxs_old]
                        # do the have overlapping atoms already?
                        if len(set(A_idxs).intersection(B_idxs)) != 0:
                            continue
                        else:
                            #they still need merging
                            # which atoms of A are closer to B center and vice versa
                            pairs = []
                            for G_idxs, ref_i in [(A_idxs, n.GetIdx()), (B_idxs, ring['atom'].GetIdx())]:
                                tm = np.take(dm, G_idxs, 0)
                                tm2 = np.take(tm, [ref_i], 1)
                                p = np.where(tm2 == np.nanmin(tm2))
                                f = G_idxs[int(p[0][0])]
                                tm2[p[0][0], :] = np.ones(tm2.shape[1]) * float('inf')
                                p = np.where(tm2 == np.nanmin(tm2))
                                s = G_idxs[int(p[1][0])]
                                pairs.append((f,s))
                            # now determine which are closer
                            if dm[pairs[0][0], pairs[1][0]] < dm[pairs[0][0], pairs[1][1]]:
                                mergituri.append((pairs[0][0], pairs[1][0]))
                                mergituri.append((pairs[0][1], pairs[1][1]))
                            else:
                                mergituri.append((pairs[0][0], pairs[1][1]))
                                mergituri.append((pairs[0][1], pairs[1][0]))
        return mergituri

    def _delete_collapsed(self, mol: Chem.RWMol):
        for a in reversed(range(mol.GetNumAtoms())):
            if mol.GetAtomWithIdx(a).GetIntProp('_ori_i') == -1:
                mol.RemoveAtom(a)

    def _mark_neighbors(self, mol, rings):
        # optional to help quinones
        for ring in rings:
            for n, bt in zip(ring['neighbors'], ring['bonds']):
                try:
                    ni = self._get_new_index(mol, n, search_collapsed=False)
                    mol.GetAtomWithIdx(ni).SetProp('_ring_bond', bt)
                except:
                    pass

    def _connenct_ring(self, mol, rings):
        # get ring members.
        old_ringers = []
        new_ringers = []
        for ring in rings:
            for i in range(len(ring['elements'])):
                new_i = self._get_new_index(mol, ring['ori_is'][i], search_collapsed=False)
                old_ringers.append(i)
                new_ringers.append(new_i)
        # keep track of new atoms
        for i in new_ringers:
            mol.GetAtomWithIdx(i).SetIntProp('expanded', 1)
        # fix ring neighbours
        for ring in rings:
            for i in range(len(ring['elements'])):
                d = self._get_expansion_for_atom(ring, i)
                new_i = self._get_new_index(mol, d['ori_i'], search_collapsed=False)
                # restore ring neighbours
                for old_neigh, bond in zip(d['neighbor'], d['bond']):
                    if old_neigh not in ring['ori_is']:
                        continue
                    bt = getattr(Chem.BondType, bond)
                    new_neigh = self._get_new_index(mol, old_neigh, search_collapsed=False)
                    present_bond = mol.GetBondBetweenAtoms(new_i, new_neigh)
                    if present_bond is None:
                        mol.AddBond(new_i, new_neigh, bt)
                    elif present_bond.GetBondType().name != bond:
                        if self._debug_draw:
                            print(
                                f'bond between {new_i} {new_neigh} exists already (has {present_bond.GetBondType().name} expected {bt})')
                        present_bond.SetBondType(bt)
                    else:
                        if self._debug_draw:
                            print(f'bond between {new_i} {new_neigh} exists already ' + \
                                  f'(has {present_bond.GetBondType().name} expected {bt})')
                        pass

    def _infer_bonding_by_proximity(self, mol):
        # fix neighbours
        # this should not happen. But just in case!
        while True:
            ringsatoms = self._get_ring_info(mol)
            for ringA, ringB in itertools.combinations(ringsatoms, 2):
                n = mol.GetNumAtoms()
                self.absorb_overclose(mol, ringA, ringB, cutoff=1.)
                if n != mol.GetNumAtoms():
                    break
            else:
                break
        new_ringers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp('expanded')]
        self.absorb_overclose(mol, new_ringers)
        new_ringers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp('expanded')]
        self.join_overclose(mol, new_ringers)
        # special case: x0749. bond between two rings
        self.join_rings(mol)

    def _get_ring_info(self, mol):
        """
        Indentical copy of fx in rectifier...
        you cannot get ring info on an unsanitized mol. Ironically I need ring info for sanitization
        :return:
        """
        mol2 = Chem.Mol(mol)
        for bond in mol2.GetBonds():
            bond.SetBondType(Chem.BondType.UNSPECIFIED)
        for atom in mol2.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetAtomicNum(0)
        Chem.SanitizeMol(mol2)
        return mol2.GetRingInfo().AtomRings()

    def join_rings(self, mol: Chem.RWMol, cutoff=1.8):
        rings = self._get_ring_info(mol)
        dm = Chem.Get3DDistanceMatrix(mol)
        for ringA, ringB in itertools.combinations(rings, 2):
            if not self._are_rings_bonded(mol, ringA, ringB):
                mini = np.take(dm, ringA, 0)
                mini = np.take(mini, ringB, 1)
                d = np.nanmin(mini)
                if d < cutoff:
                    p = np.where(mini == d)
                    f = ringA[int(p[0][0])]
                    s = ringB[int(p[1][0])]
                    mol.AddBond(f, s)



    def _are_rings_bonded(self, mol: Chem.Mol, ringA: Tuple[int], ringB: Tuple[int]):
        for i in ringA:
            for j in ringB:
                if mol.GetBondBetweenAtoms(i, j) is not None:
                    return True
        else:
            return False

    def _is_triangle(self, first: Chem.Atom, second: Chem.Atom) -> bool:
        """
        Get bool of whether two atoms share a common neighbor. Ie. joining them would make a triangle.

        :param first:
        :param second:
        :return:
        """
        get_neigh_idxs = lambda atom: [neigh.GetIdx() for neigh in atom.GetNeighbors()]
        f_neighs = get_neigh_idxs(first)
        s_neighs = get_neigh_idxs(second)
        return not set(f_neighs).isdisjoint(set(s_neighs))

    def _is_square(self, first: Chem.Atom, second: Chem.Atom) -> bool:
        """
        Get bool of whether two atoms share a common over-neighbor. Ie. joining them would make a square.

        :param first:
        :param second:
        :return:
        """
        for third in [neigh for neigh in second.GetNeighbors() if neigh.GetIdx() != first.GetIdx()]:
            if self._is_triangle(first, third) is True:
                return True
        else:
            return False

    def _is_connected_warhead(self, atom, anchor_atom):
        if not atom.HasProp('_Warhead'):
            return False
        elif atom.GetBoolProp('_Warhead') == False:
            return False
        else:
            frags = Chem.GetMolFrags(atom.GetOwningMol())
            if len(frags) == 1:
                return True
            else:
                for frag in frags:
                    if atom.GetIdx() in frag and anchor_atom.GetIdx() in frag:
                        return True
                    elif atom.GetIdx() in frag:
                        return False
                    else:
                        pass
                else:
                    raise ValueError('I do not think this is possible.')

    def join_overclose(self, mol: Chem.RWMol, to_check, cutoff=2.2): # was 1.8
        """
        Cutoff is adapted to element.

        :param mol:
        :param to_check: list of atoms indices that need joining (but not to each other)
        :param cutoff: CC bond
        :return:
        """
        pt = Chem.GetPeriodicTable()
        dm = Chem.Get3DDistanceMatrix(mol)
        for i in to_check:
            atom_i = mol.GetAtomWithIdx(i)
            for j, atom_j in enumerate(mol.GetAtoms()):
                # calculate cutoff if not C-C
                if atom_i.GetSymbol() == '*' or atom_j.GetSymbol() == '*':
                    ij_cutoff = cutoff
                elif atom_i.GetSymbol() == 'C' and atom_j.GetSymbol() == 'C':
                    ij_cutoff = cutoff
                else:
                    ij_cutoff = cutoff - 1.36 + sum([pt.GetRcovalent(atom.GetAtomicNum()) for atom in (atom_i, atom_j)])
                # determine if to join
                if i == j or j in to_check:
                    continue
                elif dm[i, j] > ij_cutoff:
                    continue
                elif self._is_triangle(atom_i, atom_j):
                    continue
                elif self._is_square(atom_i, atom_j):
                    continue
                elif self._is_connected_warhead(atom_j, atom_i):
                    continue
                else:
                    present_bond = mol.GetBondBetweenAtoms(i, j)
                    if atom_j.HasProp('_ring_bond'):
                        bt = getattr(Chem.BondType, atom_j.GetProp('_ring_bond'))
                    else:
                        bt = None
                    if present_bond is not None and bt is None:
                        pass # exists
                    elif present_bond is not None and present_bond.GetBondType() is None:
                        present_bond.SetBondType(Chem.BondType.SINGLE)
                    elif present_bond is not None and present_bond.GetBondType() is not None and bt.name == present_bond.GetBondType().name:
                        pass # exists and has correct bond
                    elif present_bond is not None:
                        present_bond.SetBondType(bt)
                    elif len(atom_i.GetNeighbors()) <= 2 and atom_i.GetIsAromatic():
                        mol.AddBond(i, j, Chem.BondType.SINGLE)
                    elif len(atom_i.GetNeighbors()) <= 3 and not atom_i.GetIsAromatic():
                        if bt is None:
                            bt = Chem.BondType.SINGLE
                        mol.AddBond(i, j, bt)
                    else:
                        pass # too bonded already!

    def absorb_overclose(self, mol, to_check=None, to_check2=None, cutoff:float=1.):
        # to_check list of indices to check and absorb into, else all atoms are tested
        dm = Chem.Get3DDistanceMatrix(mol)
        morituri = []
        if to_check is None:
            to_check = list(range(mol.GetNumAtoms()))
        if to_check2 is None:
            to_check2 = list(range(mol.GetNumAtoms()))
        for i in to_check:
            if i in morituri:
                continue
            for j in to_check2:
                if i == j or j in morituri:
                    continue
                elif dm[i, j] < cutoff:
                    self._absorb(mol, i, j)
                    morituri.append(j)
                else:
                    pass
        # kill morituri
        for i in sorted(morituri, reverse=True):
            mol.RemoveAtom(i)
        return len(morituri)

    def _absorb(self, mol, i: int, j: int):
        log.debug(f'Absorbing atom {i} with {j}')
        absorbiturum = mol.GetAtomWithIdx(j)
        for neighbor in absorbiturum.GetNeighbors():
            n = neighbor.GetIdx()
            bt = mol.GetBondBetweenAtoms(j, n).GetBondType()
            # mol.RemoveBond(j, n)
            if i == n:
                pass
            elif mol.GetBondBetweenAtoms(i, n) is None:
                mol.AddBond(i, n, bt)
            else:
                pass  # don't bother checking if differ


    def _restore_original_bonding(self, mol: Chem.RWMol, rings) -> None:
        to_be_waited_for = []
        for ring in rings:
            for i in range(len(ring['elements'])):
                d = self._get_expansion_for_atom(ring, i)
                new_i = self._get_new_index(mol, d['ori_i'], search_collapsed=False)
                for old_neigh, bond in zip(d['neighbor'], d['bond']):
                    bt = getattr(Chem.BondType, bond)
                    try:
                        new_neigh = self._get_new_index(mol, old_neigh, search_collapsed=False)
                        present_bond = mol.GetBondBetweenAtoms(new_i, new_neigh)
                        if present_bond is None:
                            mol.AddBond(new_i, new_neigh, bt)
                        elif present_bond.GetBondType().name != bond:
                            if self._debug_draw:
                                print(
                                    f'bond between {new_i} {new_neigh} exists already (has {present_bond.GetBondType().name} expected {bt})')
                            present_bond.SetBondType(bt)
                        else:
                            if self._debug_draw:
                                print(f'bond between {new_i} {new_neigh} exists already ' + \
                                      f'(has {present_bond.GetBondType().name} expected {bt})')
                            pass
                    except ValueError:
                        if self._debug_draw:
                            print(f"The neighbour {old_neigh} of {d['ori_i']} with {bt} does not yet exist")
                        to_be_waited_for.append((new_i, old_neigh, bt))
        for new_i, old_neigh, bt in to_be_waited_for:
            try:
                new_neigh = self._get_new_index(mol, old_neigh, name_restriction=mol.GetAtomWithIdx(new_i).GetProp('_ori_name'))
                if self._debug_draw:
                    print(f'{old_neigh} was missing, but has appeared since as {new_neigh}')
                if not mol.GetBondBetweenAtoms(new_i, new_neigh):
                    mol.AddBond(new_i, new_neigh, bt)
            except (KeyError, ValueError) as err:
                warn(str(err))

    # =========== Collapse =============================================================================================

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
                        else:
                            for other_center_i in old2center[neigh]:
                                if center_i != other_center_i:
                                    if not mol.GetBondBetweenAtoms(center_i, other_center_i):
                                        mol.AddBond(center_i, other_center_i, bt)
                                    break
                            else:
                                raise ValueError(f'Cannot find what {neigh} became')
        for i in sorted(set(morituri), reverse=True):
            mol.RemoveAtom(self._get_new_index(mol, i))
        return mol.GetMol()

    # =========== Misc =================================================================================================

    def _print_stored(self, mol: Chem.Mol):
        print('Idx', 'OriIdx', 'OriName')
        for atom in mol.GetAtoms():
            try:
                if atom.GetIntProp('_ori_i') == -1:
                    print(atom.GetIdx(),
                          atom.GetSymbol(),
                          atom.GetIntProp('_ori_i'),
                          atom.GetProp('_ori_name'),
                          atom.GetProp('_ori_is'),
                          atom.GetProp('_neighbors')
                          )
                else:
                    print(atom.GetIdx(),
                          atom.GetSymbol(),
                          atom.GetIntProp('_ori_i'),
                          atom.GetProp('_ori_name'))
            except KeyError:
                print(atom.GetIdx(),
                      atom.GetSymbol(), '!!!!')

    def _get_ori_i(self, mol: Chem.Mol, include_collapsed=True):
        indices = [atom.GetIntProp('_ori_i') for atom in mol.GetAtoms()]
        if include_collapsed:
            for atom in self._get_collapsed_atoms(mol):
                indices.extend(json.loads(atom.GetProp('_ori_is')))
        else:
            pass
        return indices

    def _get_collapsed_atoms(self, mol: Chem.Mol) -> List[Chem.Atom]:
        return [atom for atom in mol.GetAtoms() if atom.GetIntProp('_ori_i') == -1]

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

    def _get_mystery_ori_i(self, mol: Chem.Mol):
        """
        Collapsed rings may remember neighbors that do not exist any more...

        :param mol:
        :return:
        """
        present = self._get_ori_i(mol)
        absent = []
        for atom in self._get_collapsed_atoms(mol):
            neighss = json.loads(atom.GetProp('_neighbors'))
            absent.extend([n for neighs in neighss for n in neighs if n not in present])
        return absent

    def _is_present(self, mol, i):
        try:
            self._get_new_index(mol, i, search_collapsed=False)
            # raises value error.
            return True
        except ValueError as err:  # no atom is present (actually the default)
            return False