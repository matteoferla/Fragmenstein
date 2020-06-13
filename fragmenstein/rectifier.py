########################################################################################################################

__doc__ = \
    """
Fix issue in auto-merging
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################


from rdkit import Chem
from typing import Optional, List, Tuple
import warnings, random


# when hits are combined they can result in odd valence and other issues.


class Rectifier:
    """
    Checks whether the valence is right and corrects by either shifting the element or adding a charge.
    With the following exceptions:

    * Texas carbon -> Sulfur
    * Hydrogen -> Fluoride shifted downwards
    * Carbon in aromatic ring -> single bond ring

    :param mol:
    :param self.valence_correction: 'element' or 'charge'
    :return:
    """

    def __init__(self, mol: Chem.Mol, valence_correction: str = 'element', debug: bool = False):
        self.valence_correction = valence_correction
        if debug:
            self.dprint = print
        else:
            self.dprint = lambda *x: None
        self.mol = mol
        self._iterations_done = 0
        self.ununspecified_bonds()
        self.triage_rings()
        self.fix_issues()
        Chem.SanitizeMol(self.mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)

    # ========= Methods that circumvent the nonsanitization ============================================================

    def _get_valence_difference(self, atom: Chem.Atom) -> int:
        pt = Chem.GetPeriodicTable()
        valence = self._get_atom_valence(atom)
        maxv = max(pt.GetValenceList(atom.GetAtomicNum()))
        return valence - maxv

    def _get_atom_valence(self, atom: Chem.Atom):
        """
        Cannot get the normal way as it cannot be sanitised.

        :param atom:
        :return:
        """
        valence = 0
        for bond in atom.GetBonds():
            valence += bond.GetBondTypeAsDouble()
        return valence - atom.GetFormalCharge()

    def _has_correct_valence(self, atom: Chem.Atom):
        return self._get_valence_difference(atom) <= 0

    def _get_ring_info(self):
        """
        you cannot get ring info on an unsanitized mol. Ironically I need ring info for sanitization
        :return:
        """
        mol2 = Chem.Mol(self.mol)
        for bond in mol2.GetBonds():
            bond.SetBondType(Chem.BondType.UNSPECIFIED)
        for atom in mol2.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetAtomicNum(0)
        Chem.SanitizeMol(mol2)
        return mol2.GetRingInfo().AtomRings()

    # ========= rings ==================================================================================================

    def _is_aromatic_ring(self, ring: Tuple[int]) -> bool:
        """
        :param ring: GetRingInfo().AtomRings() entry
        :return:
        """
        for i in ring:
            atom_i = self.mol.GetAtomWithIdx(i)
            for n in atom_i.GetNeighbors():
                ni = n.GetIdx()
                if ni in ring:
                    if self.mol.GetBondBetweenAtoms(i, ni).GetBondType().name == 'AROMATIC':
                        return True
        else:
            return False

    def _get_ring_neighbors(self, ring: Tuple[int]) -> List[Tuple[int, int]]:
        """
        :param ring: GetRingInfo().AtomRings() entry
        :return: list of pairs of indices that are neighbors in the ring
        """
        rns = []
        for i in ring:
            atom = self.mol.GetAtomWithIdx(i)
            for n in atom.GetNeighbors():
                ni = n.GetIdx()
                if ni in ring:
                    rns.append((i, ni))
        return rns

    def triage_rings(self):
        # upgrade to aromatic if aromatic.
        rings = self._get_ring_info()
        for ring in rings:
            if self._is_aromatic_ring(ring):
                for i, ni in self._get_ring_neighbors(ring):
                    self.mol.GetBondBetweenAtoms(i, ni).SetBondType(Chem.BondType.AROMATIC)
        # downgrade to single if non ring aromatic
        for i, atom in enumerate(self.mol.GetAtoms()):
            if any([i in r for r in rings]):
                continue  # ring
            else:  # non-ring
                for bond in atom.GetBonds():
                    if bond.GetBondType().name == 'AROMATIC':
                        self.dprint(f'donwgrading bond {i}')
                        bond.SetBondType(Chem.BondType.SINGLE)
        # aromatics
        for ring in sorted(rings, key=self._is_aromatic_ring):
            if self._is_aromatic_ring(ring): # the nonaromatic rings will be done first.
                for i in ring:
                    self.mol.GetAtomWithIdx(i).SetIsAromatic(True)
            else:
                for i in ring:
                    self.mol.GetAtomWithIdx(i).SetIsAromatic(False)


    def _get_aroma(self, atom, this_bond):
        return [b for b in atom.GetBonds() if b.GetIdx() != this_bond and b.GetBondType().name == 'AROMATIC']

    def _get_other(self, bond, these_atoms):
        others = [a for a in (bond.GetBeginAtom(), bond.GetEndAtom()) if a.GetIdx() not in these_atoms]
        if others:
            other = others[0]
            other.SetIsAromatic(False)
            return other

    def downgrade_ring(self, atom: Chem.Atom):
        ## very crappy way of doing this
        self.dprint(f'downgrading whole ring')
        atom.SetIsAromatic(False)
        for bond in atom.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
            other = self._get_other(bond, [atom.GetIdx()])
            aro = self._get_aroma(other, bond.GetIdx())
            if aro:
                aro[0].SetBondType(Chem.BondType.DOUBLE)
                doubleother = self._get_other(aro[0], [atom.GetIdx(), other.GetIdx()])
                for b in doubleother.GetBonds():
                    if b.GetBondType() == Chem.BondType.AROMATIC:
                        b.SetBondType(Chem.BondType.SINGLE)
                    neigh = self._get_other(b, [doubleother.GetIdx()])
                    if neigh:
                        neigh.SetIsAromatic(False)

    # ========= Sanitization based fixes ===============================================================================

    def fix_issues(self):
        if self._iterations_done > 10:
            return None
        problems = Chem.DetectChemistryProblems(self.mol)
        if len(problems) == 0:
            return None
        else:
            self._iterations_done += 1
        for p in problems:
            ############################################################
            if p.GetType() == 'KekulizeException':
                # plural GetAtomIndices. AtomKekulizeException, singular GetAtomIdx
                print(p.Message())
                N = self._get_nitrogens(p.GetAtomIndices())
                if len(N) > 0:
                    random.shuffle(N)  # just in case.
                    self.mol.GetAtomWithIdx(N[0]).SetNumExplicitHs(1)
            ############################################################
            elif p.GetType() == 'AtomKekulizeException' and 'non-ring atom' in p.Message():
                atom = self.mol.GetAtomWithIdx(p.GetAtomIdx())
                atom.SetIsAromatic(False)
                for bond in atom.GetBonds():
                    bond.SetBondType(Chem.BondType.SINGLE)
            ############################################################
            elif p.GetType() == 'AtomValenceException':
                i = p.GetAtomIdx()
                self.fix_valence(i)
            else:
                self.dprint('???', p.GetType(), p.Message())
        self.fix_issues()

    # ========= other helpers ==========================================================================================

    def _get_nitrogens(self, indices):
        return [i for i in indices if self.mol.GetAtomWithIdx(i).GetSymbol() == 'N']

    def ununspecified_bonds(self):
        for bond in self.mol.GetBonds():
            if bond.GetBondType().name == 'UNSPECIFIED':
                bond.SetBondType(Chem.BondType.SINGLE)

    # ========= shift/charge ===========================================================================================

    def _adjust_for_fix_valence(self, atom):
        df = self._get_valence_difference(atom)
        if self.valence_correction == 'charge':
            atom.SetFormalCharge(df)
        elif self.valence_correction == 'element':
            if atom.GetSymbol == 'H':
                atom.SetAtomicNum(8)
            n = atom.GetAtomicNum()
            if len(atom.GetNeighbors()) > 4:
                atom.SetAtomicNum(16)
            elif n - df < 6:  # C -> B no!
                for bond in atom.GetBonds():
                    bond.SetBondType(Chem.BondType.SINGLE)
            else:  # N, O, F etc.
                atom.SetAtomicNum(int(n - df))
        else:
            raise ValueError(f'self.valence_correction can only be "element"/"charge" not {self.valence_correction}.')

    def fix_valence(self, i):
        atom = self.mol.GetAtomWithIdx(i)
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        self.dprint(f'{i} {atom.GetSymbol()}: {len(atom.GetNeighbors())} bonds {self._get_atom_valence(atom)}')
        if self._has_correct_valence(atom):
            self.dprint('\tValence seems correct')
            return None
        elif atom.GetSymbol() == 'C' and len(atom.GetNeighbors()) > 4:
            self.dprint('\ttexas carbon --> S')
            atom.SetAtomicNum(16)
        elif atom.GetSymbol() == 'C' and atom.GetIsAromatic() and len(atom.GetNeighbors()) == 4:
            self.dprint('\tDowngrading ring')
            self.downgrade_ring(atom)
        elif atom.GetSymbol() == 'C':
            for bond in atom.GetBonds():
                bond.SetBondType(Chem.BondType.SINGLE)
        else:
            self._adjust_for_fix_valence(atom)
        # did it work?
        if self._has_correct_valence(atom):
            return self.mol
        else:
            return self.fix_valence(i)