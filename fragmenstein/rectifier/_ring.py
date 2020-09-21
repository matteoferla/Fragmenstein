########################################################################################################################

__doc__ = \
    """
    This add the ring fixing functionality. Fusing etc.
    Formerly part of collapse_ring.py
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
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Optional, Dict, List, Any, Tuple, Union
import numpy as np
from collections import Counter, defaultdict

from ._base import _RectifierBaseMixin

# ======================================================================================================================


class _RectifierRingMixin(_RectifierBaseMixin):

    def fix_rings(self):
        self._prevent_conjoined_ring()
        self._prevent_weird_rings()

    def _prevent_conjoined_ring(self) -> None:
        """
        This kills bridging bonds with not atoms in the bridge within rings.
        So it is bridged, fused and spiro safe.
        It removes only one bond, so andamantane/norbornane are safe.
        """
        c = Counter([i for ring in self._get_ring_info() for i in ring])
        nested = [k for k in c if c[k] >= 3]
        pairs = [(idx_a, idx_b) for idx_a, idx_b in itertools.combinations(nested, r=2) if
                 self.rwmol.GetBondBetweenAtoms(idx_a, idx_b) is not None]
        rank = sorted(pairs, key=lambda x: c[x[0]] + c[x[1]], reverse=True)
        if len(rank) > 0:
            idx_a, idx_b = rank[0]
            self.rwmol.RemoveBond(idx_a, idx_b)  # SetBoolProp('_IsRingBond') is not important
            self.journal.info(f'Zero-atom bridged ring issue: bond between {idx_a}-{idx_b} removed')
            # re-run:
            self._prevent_conjoined_ring()
        self.modifications.append(self.mol)

    def _prevent_weird_rings(self):
        ringatoms = self._get_ring_info()  # GetRingInfo().AtomRings()
        for ring_A, ring_B in itertools.combinations(ringatoms, r=2):
            shared = set(ring_A).intersection(set(ring_B))
            if len(shared) == 0:
                self.journal.debug('This molecule has some separate rings')
                pass  # separate rings
            elif len(shared) < self.atoms_in_bridge_cutoff and \
                    self.atoms_in_bridge_cutoff >= 2 \
                    and len(ring_A) == len(ring_B):
                # adamantene/norbornane/tropinone kind of thing
                self.journal.warning('This molecule has a bridge: leaving')
                pass  # ideally check if planar...
            elif len(shared) == 1:
                self.journal.debug('This molecule has a spiro bicycle')
                pass  # spiro ring.
            elif len(shared) == 2:
                self.journal.debug('This molecule has a fused ring')
                if self.rwmol.GetBondBetweenAtoms(*shared) is not None:
                    pass  # indole/naphtalene
                    small, big = sorted([ring_A, ring_B], key=lambda ring: len(ring))
                    if len(small) == 4:
                        self.journal.warning('This molecule has a benzo-azetine–kind-of-thing: expanding to indole')
                        # Chem.MolFromSmiles('C12CCCCC1CC2')
                        # benzo-azetine is likely an error: add and extra atom
                        a, b = set(small).difference(big)
                        self._place_between(a, b)
                    elif len(small) == 3:
                        self.journal.warning('This molecule has a benzo-cyclopropane–kind-of-thing: expanding to indole')
                        # Chem.MolFromSmiles('C12CCCCC1C2')
                        # benzo-cyclopronane is actually impossible at this stage.
                        a = list(set(small).difference(big))[0]
                        for b in shared:
                            self._place_between(a, b)
                    else:
                        pass  # indole and nathalene
                elif (len(ring_A), len(ring_B)) == (6, 6):
                    raise Exception('This is utterly impossible')
                else:
                    self.journal.warning(f'mysterious ring system {len(ring_A)} + {len(ring_B)}')
                    pass  # ????
            elif len(shared) < self.atoms_in_bridge_cutoff:
                # adamantene/norbornane/tropinone kind of thing
                self.journal.warning('This molecule has a bridge: leaving')
                pass  # ideally check if planar...
            else:
                self.journal.warning('This molecule has a bridge that will be removed')
                self._prevent_bridge_ring(ring_A)
                # start from scratch.
                self._prevent_weird_rings()
        self.modifications.append(self.mol)

    # ===== Dep of prevent weird rings =================================================================================

    def _place_between(self, a: int, b: int, aromatic: Optional[bool] = None, atomic_number: int = 6) -> None:
        """
        Places an C atom, possibly of type aromatic, between atom of index a, and of b.

        :param a: index of atom A
        :param b: index of atom B
        :param aromatic: bool of aromaticity (False = Single, None = copy, True = aromatic)
        :param atomic_number: Carbon is 6.
        :return:
        """
        oribond = self.rwmol.GetBondBetweenAtoms(a, b)
        if oribond is None:
            self.journal.critical(f'FAIL. There should be a bond btween {a} and {b}')
            return None  # fail
        elif aromatic is True:
            bt = Chem.BondType.AROMATIC
        elif aromatic is False:
            bt = Chem.BondType.SINGLE
        else:
            bt = oribond.GetBondType()
        idx = self.rwmol.AddAtom(Chem.Atom(atomic_number))
        neoatom = self.rwmol.GetAtomWithIdx(idx)
        atom_a = self.rwmol.GetAtomWithIdx(a)
        atom_b = self.rwmol.GetAtomWithIdx(b)
        if aromatic:
            neoatom.SetIsAromatic(True)
            atom_a.SetIsAromatic(True)
            atom_b.SetIsAromatic(True)
        # prevent constraints
        neoatom.SetBoolProp('_Novel', True)
        atom_a.SetBoolProp('_Novel', True)
        atom_b.SetBoolProp('_Novel', True)
        # fix position
        conf = self.rwmol.GetConformer()
        pos_A = conf.GetAtomPosition(a)
        pos_B = conf.GetAtomPosition(b)
        x = pos_A.x / 2 + pos_B.x / 2
        y = pos_A.y / 2 + pos_B.y / 2
        z = pos_A.z / 2 + pos_B.z / 2
        conf.SetAtomPosition(idx, Point3D(x, y, z))
        # fix bonds
        self.rwmol.RemoveBond(a, b)
        self.rwmol.AddBond(a, idx, bt)
        self.rwmol.AddBond(b, idx, bt)

    def _prevent_bridge_ring(self, examplar: Tuple[int]) -> None:
        # examplar is ring
        ringatoms = self._get_ring_info()  # GetRingInfo().AtomRings()
        ringatoms = [ring for ring in ringatoms if set(ring).intersection(examplar)]
        ring_idx = list(range(len(ringatoms)))
        shared_count = {}
        for ra, rb in itertools.combinations(ring_idx, r=2):
            shared_count[(ra, rb)] = len(set(ringatoms[ra]).intersection(set(ringatoms[rb])))
        if len(shared_count) == 0:
            return None
        ra, rb = list(shared_count.keys())[0]
        shared = list(set(ringatoms[ra]).intersection(ringatoms[rb]))
        has_bond = lambda a, b: self.rwmol.GetBondBetweenAtoms(a, b) is not None
        pairs = [(a, b) for a, b in itertools.combinations(shared, r=2) if has_bond(a, b)]
        c = Counter([i for pair in pairs for i in pair])
        ring_A, ring_B = ringatoms[ra], ringatoms[rb]
        small, big = sorted([ring_A, ring_B], key=lambda ring: len(ring))
        inners = [i for i in c if c[i] > 1]
        x = list(set(shared).difference(inners))
        if len(x) != 2:
            self.journal.critical(
                f'This is impossible. {ringatoms} share {shared} with {inners} in the inside and {x} on the edge?')
            return None
        a, b = x
        if len(big) > 6:
            self.journal.warning(f'Removing {len(inners)} bridging atoms and replacing with fused ring')
            # bond the vertices
            bt = Chem.BondType.SINGLE  # ???
            if self.rwmol.GetBondBetweenAtoms(a, b) is None:
                self.rwmol.AddBond(a, b, bt)
            else:
                self.journal.warning('This is really odd! Why is there a bond already??')
            # remove the middle atoms.
            for i in sorted(inners, reverse=True):
                self.rwmol.RemoveAtom(i)
        else:
            self.journal.warning(f'Shriking the smaller ring to change from bridged to fused.')
            # get the neighbour in the small atom to a vertex.
            neighs = [neigh for neigh in self.rwmol.GetAtomWithIdx(a).GetNeighbors() if
                      neigh.GetIdx() not in shared and neigh.GetIdx() in small]
            neigh = sorted(neighs, key=lambda atom: atom.GetSymbol() != 'C')[0]
            bt = self.rwmol.GetBondBetweenAtoms(a, neigh.GetIdx()).GetBondType()
            self.rwmol.RemoveBond(a, neigh.GetIdx())
            new_neigh = [neigh for neigh in self.rwmol.GetAtomWithIdx(a).GetNeighbors() if neigh.GetIdx() in shared][0]
            self.rwmol.AddBond(neigh.GetIdx(), new_neigh.GetIdx(), bt)
            neigh.SetBoolProp('_Novel', True)
            new_neigh.SetBoolProp('_Novel', True)
            self.rwmol.GetAtomWithIdx(a).SetBoolProp('_Novel', True)


