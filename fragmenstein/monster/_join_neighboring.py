########################################################################################################################
__doc__ = \
    """
This is inherited by both place and combine via _MonsterMerge
    """

########################################################################################################################

from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D
from typing import Tuple, List, Dict, Optional, Union
import numpy as np
from warnings import warn
from .bond_provenance import BondProvenance
from ._communal import _MonsterCommunal

class _MonsterJoinNeigh(_MonsterCommunal):
    def join_neighboring_mols(self, mol_A: Chem.Mol, mol_B: Chem.Mol):
        """
        Joins two molecules by first calling _find_closest to find closest.
        That method does all the thinking.
        then by calling _join_atoms.

        :param mol_A:
        :param mol_B:
        :return:
        """
        # get closets atoms
        combo, candidates = self._find_all_closest(mol_A, mol_B)  # _find_all_closest is in communal
        anchor_A, anchor_B, distance = candidates[0]
        mol = self._join_atoms(combo, anchor_A, anchor_B, distance, linking=True)
        for anchor_A, anchor_B, distance in candidates[1:]:
            mol = self._join_atoms(combo, anchor_A, anchor_B, distance, linking=False)
        mol.SetProp('_Name', mol_A.GetProp('_Name') + '~' + mol_B.GetProp('_Name'))
        return mol


    def _join_atoms(self,
                    combo: Chem.RWMol,
                    anchor_A: int,
                    anchor_B: int,
                    distance: float,
                    linking: bool=True):
        """
        extrapolate positions between. by adding linkers if needed.
        """
        conf = combo.GetConformer()
        pos_A = conf.GetAtomPosition(anchor_A)
        pos_B = conf.GetAtomPosition(anchor_B)
        n_new = int(round(distance / 1.22) - 1)
        xs = np.linspace(pos_A.x, pos_B.x, n_new + 2)[1:-1]
        ys = np.linspace(pos_A.y, pos_B.y, n_new + 2)[1:-1]
        zs = np.linspace(pos_A.z, pos_B.z, n_new + 2)[1:-1]

        # correcting for ring marker atoms
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
            self.journal.debug(f'Two ring atoms detected to be close. Joining for now.'+
                               ' They will be bonded/fused/spiro afterwards')
        # check if valid.
        if distance > self.joining_cutoff:
            msg = f'Atoms {anchor_A}+{anchor_B} are {distance} Å away. Cutoff is {self.joining_cutoff}.'
            self.journal.warning(msg)
            raise ConnectionError(msg)
        # place new atoms
        self.journal.debug(f'Molecules will be joined via atoms {anchor_A}+{anchor_B} ({distance} Å) via the addition of {n_new} atoms.')
        previous = anchor_A
        if linking is False and n_new > 0:
            self.journal.warning(f'Was going to bond {anchor_A} and {anchor_B} but reconsidered.')
        elif linking is True and n_new <= 0:
            combo.AddBond(previous, anchor_B, Chem.BondType.SINGLE)
            new_bond = combo.GetBondBetweenAtoms(previous, anchor_B)
            BondProvenance.set_bond(new_bond, 'main_novel')
        elif linking is False and n_new <= 0:
            combo.AddBond(previous, anchor_B, Chem.BondType.SINGLE)
            new_bond = combo.GetBondBetweenAtoms(previous, anchor_B)
            BondProvenance.set_bond(new_bond, 'other_novel')
        elif linking is True and n_new > 0:
            for i in range(n_new):
                # make oxygen the first and last bridging atom.
                if i == 0 and combo.GetAtomWithIdx(anchor_A).GetSymbol() == 'C':
                    new_atomic = 8
                elif i > 2 and i == n_new -1 and combo.GetAtomWithIdx(anchor_B).GetSymbol() == 'C':
                    new_atomic = 8
                else:
                    new_atomic = 6
                idx = combo.AddAtom(Chem.Atom(new_atomic))
                new = combo.GetAtomWithIdx(idx)
                new.SetBoolProp('_Novel', True)
                new.SetIntProp('_ori_i', 999)
                conf.SetAtomPosition(idx, Point3D(float(xs[i]), float(ys[i]), float(zs[i])))
                combo.AddBond(idx, previous, Chem.BondType.SINGLE)
                new_bond = combo.GetBondBetweenAtoms(idx, previous)
                BondProvenance.set_bond(new_bond, 'linker')
                previous = idx
            combo.AddBond(previous, anchor_B, Chem.BondType.SINGLE)
            new_bond = combo.GetBondBetweenAtoms(previous, anchor_B)
            BondProvenance.set_bond(new_bond, 'linker')
        else:
            raise ValueError('Impossible')
        return combo.GetMol()

    def get_largest_fragment(self, mol):
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        frags = sorted(frags, key=lambda mol: mol.GetNumAtoms(), reverse=True)
        discarded_origins = {a.GetProp('_ori_name') for m in frags[1:] for a in m.GetAtoms() if a.HasProp('_ori_name')}
        self.unmatched.extend(discarded_origins)
        return frags[0]
