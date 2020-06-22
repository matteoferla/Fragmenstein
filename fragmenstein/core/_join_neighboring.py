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
        Joins two molecules by first calling _find_closest to find closest
        then by calling _join_atoms

        :param mol_A:
        :param mol_B:
        :return:
        """
        # offset to avoid clashes in bond by history mode
        self._offset_collapsed_ring(mol_B)
        self._offset_origins(mol_B)
        # get closets atoms
        combo, anchor_A, anchor_B, d = self._find_closest(mol_A, mol_B)
        mol = self._join_atoms(combo, anchor_A, anchor_B, d)
        mol.SetProp('_Name', mol_A.GetProp('_Name') + '~' + mol_B.GetProp('_Name'))
        return mol

    def _find_closest(self, mol_A: Chem.Mol, mol_B: Chem.Mol) -> Tuple[Chem.RWMol, int, int, float]:
        """
        first step in join_neighboring_mols

        :param mol_A:
        :param mol_B:
        :return:
        """
        combo = Chem.RWMol(Chem.CombineMols(mol_A, mol_B))
        dm = Chem.Get3DDistanceMatrix(combo)
        A_idxs = np.arange(mol_A.GetNumAtoms())
        B_idxs = np.arange(mol_B.GetNumAtoms()) + mol_A.GetNumAtoms()
        tm = np.take(dm, A_idxs, 0)
        tm2 = np.take(tm, B_idxs, 1)
        d = float('inf')

        def is_fullbonded(atom):
            if atom.GetIsAromatic() and len(atom.GetNeighbors()) > 2:
                return True
            elif len(atom.GetNeighbors()) > 3:
                return True
            else:
                return False

        def is_ring_atom(atom):
            if atom.GetIsAromatic():
                return True
            elif atom.HasProp('_ori_i') and atom.GetIntProp('_ori_i') == -1:
                if atom.HasProp('_bonds') and 'AROMATIC' in atom.GetProp('_bonds'):
                    return True  # technically it could be non-aromatic (ring fusion).
                else:
                    return False
            else:
                return False

        def is_warhead_marked(atom):
            return atom.HasProp('_Warhead') and atom.GetBoolProp('_Warhead') == True

        previous = None
        while d == float('inf'):
            d = np.nanmin(tm2)
            p = np.where(tm2 == d)
            anchor_A = int(A_idxs[p[0][0]])
            anchor_B = int(B_idxs[p[1][0]])
            atom_A = combo.GetAtomWithIdx(anchor_A)
            atom_B = combo.GetAtomWithIdx(anchor_B)
            def outstrike_A(p):
                tm2[p[0][0], :] = np.ones(tm2.shape[1]) * float('inf')
            def outstrike_B(p):
                tm2[:, p[1][0]] = np.ones(tm2.shape[0]) * float('inf')

            if d == float('inf') and previous is not None:
                log.debug(f'run out of options for distance bonding, using previous')
                d, anchor_A, anchor_B = previous
            elif d == float('inf'):
                raise ConnectionError('This is impossible. Previous is absent??')
            elif is_warhead_marked(atom_A):
                log.debug(f'Atom A ({anchor_A}) is warhead.')
                # previous = (d, anchor_A, anchor_B) # forbid.
                d = float('inf')
                outstrike_A(p)
            elif is_warhead_marked(atom_B):
                log.debug(f'Atom B ({anchor_B}) is warhead.')
                # previous = (d, anchor_A, anchor_B) # forbid.
                d = float('inf')
                outstrike_B(p)
            elif is_fullbonded(atom_A):
                log.debug(f'Atom A ({anchor_A}) is already at full bond allowance')
                previous = (d, anchor_A, anchor_B)
                d = float('inf')
                outstrike_A(p)
            elif is_fullbonded(atom_B):
                log.debug(f'Atom B ({anchor_B}) is already at full bond allowance')
                previous = (d, anchor_A, anchor_B)
                d = float('inf')
                outstrike_B(p)
            elif is_ring_atom(atom_A) and previous is None:
                log.debug(f'Atom A ({anchor_A}) is a ring. Don\'t really want to connect to that')
                previous = (d, anchor_A, anchor_B)
                d = float('inf')
                outstrike_A(p)
            elif is_ring_atom(atom_B) and previous is None:
                log.debug(f'Atom B ({anchor_B}) is a ring. Don\'t really want to connect to that')
                previous = (d, anchor_A, anchor_B)
                d = float('inf')
                outstrike_B(p)
            elif previous is not None and previous[0] - d > 0.5:
                log.info(f'going with the ring then, the next one is too far.')
                d, anchor_A, anchor_B = previous
            else:
                if self._debug_draw:
                    print(is_fullbonded(atom_A), is_fullbonded(atom_B), is_ring_atom(atom_A), is_ring_atom(atom_B),
                          previous)
                pass
        return combo, anchor_A, anchor_B, d

    def _join_atoms(self, combo: Chem.RWMol, anchor_A: int, anchor_B: int, d: float):
        # extrapolate positions between
        conf = combo.GetConformer()
        pos_A = conf.GetAtomPosition(anchor_A)
        pos_B = conf.GetAtomPosition(anchor_B)
        n_new = int(round(d / 1.22) - 1)
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
            d -= 1.35 + 0.2 # Arbitrary + 0.2 to compensate for the ring not reaching (out of plane).
            n_new -= 1
            xs = xs[1:]
            ys = ys[1:]
            zs = zs[1:]

        if is_ring_atom(anchor_B):
            d -= 1.35 + 0.2 # Arbitrary + 0.2 to compensate for the ring not reaching  (out of plane).
            n_new -= 1
            xs = xs[:-1]
            ys = ys[:-1]
            zs = zs[:-1]

        # notify that things could be leary.
        if d < 0:
            log.info(f'Two ring atoms detected to be too close. Joining for now. They should be merged.')
        # place new atoms
        if d < self.joining_cutoff:
            log.debug(f'Molecules will be joined via atoms {anchor_A}+{anchor_B} ({d} Å) via the addition of {n_new} atoms.')
            if self._debug_draw:
                print(f'Adding {n_new} atoms between {anchor_A} and {anchor_B} ({d} jump)')
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
            msg = f'Atoms {anchor_A}+{anchor_B} are {d} Å away. Cutoff is {self.joining_cutoff}.'
            log.warning(msg)
            raise ConnectionError(msg)
        return combo.GetMol()