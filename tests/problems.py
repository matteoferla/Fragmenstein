import unittest

# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem

from fragmenstein import Monster, Igor
from fragmenstein.mpro import MProVictor
from typing import Union
# ======================================================================================================================
class UnresolvedProblems(unittest.TestCase):

    def setUp(self):
        Igor.init_pyrosetta()

    def assertMolEqual(self, mol1: Union[str, Chem.Mol], mol2: Union[str, Chem.Mol]):
        """
        For now I am comparing SMILES (tut tut), but I will fix it one day...
        """
        smiles1: str = mol1 if isinstance(mol1, str) else Chem.MolToSmiles(AllChem.RemoveAllHs(mol1))
        smiles2: str = mol1 if isinstance(mol2, str) else Chem.MolToSmiles(AllChem.RemoveAllHs(mol2))
        self.assertEqual(smiles1, smiles2)

    def test_recto_fail_A(self):
        """This used to fail."""
        MProVictor.monster_throw_on_discard = True
        victor = MProVictor.from_hit_codes(hit_codes=['x11612', 'x11475'])
        victor.combine()
        self.assertEqual(victor.error_msg, '', victor.error_msg)

    def test_supplementary1_to_recto_fail_A(self):
        """
        This was meant to test the above, but it works fine.
        :return:
        """
        # make hits
        chlorotoluene = Chem.MolFromSmiles('c1c(Cl)c(C)ccc1')
        AllChem.EmbedMolecule(chlorotoluene)
        chlorotoluene.SetProp('_Name', 'orthochlorotoluene')
        toluene = Chem.RWMol(chlorotoluene)
        toluene.RemoveAtom(2)
        Chem.SanitizeMol(toluene)
        toluene.SetProp('_Name', 'toluene')
        chlorobutane = Chem.RWMol(chlorotoluene)
        for n in range(chlorobutane.GetNumAtoms() - 1, 4, -1):
            chlorobutane.RemoveAtom(n)
        for atom in chlorobutane.GetAtoms():
            atom.SetIsAromatic(False)
        for bond in chlorobutane.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
        Chem.SanitizeMol(chlorobutane)
        chlorobutane.SetProp('_Name', '2-chlorobutane')
        # merge
        monster = Monster(hits=[toluene, chlorobutane]).combine(keep_all=False)
        # ======
        self.assertMolEqual(chlorotoluene, monster.positioned_mol)

    def test_supplementary2_to_recto_fail_A(self):
        """
        This was meant to test as above.
        It merges xylene with chloropentane
        :return:
        """
        #
        chloroxylene = Chem.MolFromSmiles('Cc1c(Cl)c(C)ccc1')
        AllChem.EmbedMolecule(chloroxylene)
        chloroxylene.SetProp('_Name', 'chloroxylene')
        #
        xylene = Chem.RWMol(chloroxylene)
        xylene.RemoveAtom(3)
        Chem.SanitizeMol(xylene)
        xylene.SetProp('_Name', 'xylene')
        #
        chloropentane = Chem.RWMol(chloroxylene)
        for n in range(chloropentane.GetNumAtoms() - 1, 5, -1):
            chloropentane.RemoveAtom(n)
        for atom in chloropentane.GetAtoms():
            atom.SetIsAromatic(False)
        for bond in chloropentane.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
        Chem.SanitizeMol(chloropentane)
        chloropentane.SetProp('_Name', '2-chloropentane')
        #
        monster = Monster(hits=[xylene, chloropentane]).combine(keep_all=False)
        # ======
        self.assertEqual(Chem.MolToSmiles(chloroxylene), Chem.MolToSmiles(AllChem.RemoveAllHs(monster.positioned_mol)))

    def test_longer_link(self):
        """
        This is a know failure of the algorithm.
        It does not map "properly" the case where a link is lengthened. This is problematic in general
        and MCS based approach would not work even if the bondcomparison took into account whether atom is part
        of a linear stretch.
        """
        benzylbenzene = Chem.MolFromSmiles('c1ccccc1Cc1ccccc1')
        benzylbenzene.SetProp('_Name', 'benzylbenzene')
        assert AllChem.EmbedMolecule(benzylbenzene) == 0
        # extra carbon in bridge... dibenzyl
        monster = Monster([benzylbenzene, ]).place_smiles('c1ccccc1CCc1ccccc1')
        # this _seems_ correct, but it is not as the second ring is not mapped correctly.
        self.assertEqual(benzylbenzene.GetNumAtoms(),
                         len(monster.convert_origins_to_custom_map()['benzylbenzene']),
                         'not all atoms are mapped')

if __name__ == '__main__':
    unittest.main()
