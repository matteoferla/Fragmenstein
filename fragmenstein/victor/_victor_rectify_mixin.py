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


from ._victor_base_mixin import _VictorBaseMixin

from rdkit import Chem
from typing import Optional
import warnings

# when hits are combined they can result in odd valence and other issues.

class _VictorRectifyMixin(_VictorBaseMixin):
    def rectify(self, mol:Chem.Mol, valence_correction: str = 'element'):
        """
        Checks whether the valence is right and corrects by either shifting the element or adding a charge.
        With the following exceptions:

        * Texas carbon -> Sulfur
        * Hydrogen -> Chlorine shifted downwards
        * Carbon in aromatic ring -> single bond ring

        :param mol:
        :param valence_correction: 'element' or 'charge'
        :return:
        """
        pt = Chem.GetPeriodicTable()
        for atom in mol.GetAtoms():
            if self._has_correct_valence(atom):
                continue
            elif {bond.GetBondType().name for bond in atom.GetBonds()} == 'AROMATIC':
                self.journal.debug(f'{self.long_name} - Likely false alarm {atom.GetIdx()} {atom.GetSymbol()} 3 aromatic bond.')
                continue
            elif atom.GetSymbol() == 'C' and atom.GetIsAromatic():
                self.journal.warning(f'{self.long_name} - Aromatic carbon valence issue {atom.GetIdx()} {atom.GetSymbol()}')
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.AROMATIC:
                        self._downgrade_ring(mol, atom)
            elif atom.GetSymbol() == 'C':
                # texas carbon? No sulfur. Ideally collapse of closest atom pair.
                self.journal.warning(f'{self.long_name} - Texas carbon changed to sulfur')
                atom.SetAtomicNum(16)
            else:
                self.journal.warning(f'{self.long_name} - Other valence issue {atom.GetIdx()} {atom.GetSymbol()}')
                if valence_correction == 'charge':
                    atom.SetFormalCharge(self._get_valence_difference(atom))
                elif valence_correction == 'element':
                    if atom.GetSymbol == 'H':
                        atom.SetAtomicNum(8)
                    n = atom.GetAtomicNum()
                    atom.SetAtomicNum(n - self._get_valence_difference(atom))
                else:
                    raise ValueError(f'valence_correction can only be "element" or "charge" not {valence_correction}.')
        Chem.SanitizeMol(mol, catchErrors=True)
        return mol

    def _downgrade_ring(self, mol:Chem.Mol, atom: Chem.Atom):

        def get_aroma(atom, this_bond):
            return [b for b in atom.GetBonds() if b.GetIdx() != this_bond and b.GetBondType().name == 'AROMATIC']

        def get_other(bond, these_atoms):
            others = [a for a in (bond.GetBeginAtom(), bond.GetEndAtom()) if a.GetIdx() not in these_atoms]
            if others:
                other = others[0]
                other.SetIsAromatic(False)
                return other

        # mol.GetRingInfo().AtomRings() does not work with unsanitised
        ## very crappy way of doing this
        atom.SetIsAromatic(False)
        for bond in atom.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
            other = get_other(bond, [atom.GetIdx()])
            aro = get_aroma(other, bond.GetIdx())
            if aro:
                aro[0].SetBondType(Chem.BondType.DOUBLE)
                doubleother = get_other(aro[0], [atom.GetIdx(), other.GetIdx()])
                for b in doubleother.GetBonds():
                    if b.GetBondType() == Chem.BondType.AROMATIC:
                        b.SetBondType(Chem.BondType.SINGLE)
                    neigh = get_other(b, [doubleother.GetIdx()])
                    if neigh:
                        neigh.SetIsAromatic(False)
                        
    def _get_atom_valence(self, atom: Chem.Atom):
        # cannot get the normal way or errors are raised.
        valence = 0
        for bond in atom.GetBonds():
            valence += bond.GetBondTypeAsDouble()
        return valence - atom.GetFormalCharge()
    
    def _get_valence_difference(self, atom: Chem.Atom):
        pt = Chem.GetPeriodicTable()
        valence = self._get_atom_valence(atom)
        maxv = max(pt.GetValenceList(atom.GetAtomicNum()))
        return valence - maxv

    def _has_correct_valence(self, atom: Chem.Atom):
        return self._get_valence_difference(atom) <= 0