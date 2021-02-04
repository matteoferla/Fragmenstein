import unittest, os
# ======================================================================================================================
import pyrosetta

pyrosetta.init(
    extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')
# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem

from fragmenstein import Monster, Victor, Igor, Rectifier
from fragmenstein.mpro import MProVictor
from typing import *


# ======================================================================================================================


class MProTargetTester(unittest.TestCase):

    def test_easy(self):
        """
        To a **human** this looks easy. x0692 is a red herring and the other two make the molecule.
        As currently written this will fail.

        :return:
        """
        # PAU-UNI-52c0427f-1
        MProVictor.quick_renanimation = True
        victor = MProVictor.from_hit_codes(smiles='CCNc1ncc(C#N)cc1CN1CCN(C(=O)C*)CC1',
                                           hit_codes=['x0692', 'x0305', 'x1249'],
                                           long_name='2_ACL')
        self.assertEqual(victor.error, '', victor.error)
        self.assertIsNotNone(victor.minimised_mol, 'Failed minimisation')
        msg = f'x1249 is the red herring, prediction: {victor.monster.unmatched}, ' + \
              f'while x0305 and x0692 the true inspirations {victor.monster.matched}'
        self.assertIn('x1249', victor.monster.unmatched, msg)
        self.assertIn('x0305', victor.monster.matched, msg)

        victor.make_pse()

    def test_nasty(self):
        """
        The human suggested a lot of novel groups.
        'x0540' is a really odd inspiration. Three atoms are conserved. the rest aren't.

        :return:
        """
        MProVictor.quick_renanimation = True
        victor = MProVictor.from_hit_codes(smiles='*CCC(=O)N1CC(CCN(C(=O)Nc2c(C)ncc(C)c2CCN2CCOCC2)c2cc(C)ccn2)C1',
                                           hit_codes=['x0434', 'x0540'],
                                           long_name='AGN-NEW-5f0-1_ACR1')
        self.assertEqual(str(victor.error), '', str(victor.error))
        self.assertIsNotNone(victor.minimised_mol, 'Failed minimisation')
        self.assertEqual(len(victor.monster.unmatched), 0,
                         f'Both were correct but {victor.monster.unmatched} was discarded')
        victor.make_pse()

    def test_incorrect(self):
        """
        This case has two hits that are not ideal. One, x0995, totally barney.
        These will be rejected.

        :return:
        """
        MProVictor.quick_renanimation = True
        victor = MProVictor.from_hit_codes(smiles='*C(=N)CN1CCN(Cc2ccc(-c3cc(CC)ncn3)c(F)c2)CC1',
                                           hit_codes='x0692,x0770,x0995'.split(','),
                                           long_name='BEN-VAN-c98-4')
        self.assertEqual(victor.error, '', victor.error)
        self.assertIsNotNone(victor.minimised_mol, 'Failed minimisation')
        victor.make_pse()
        self.assertIn('x0995', victor.monster.unmatched)

    def test_pentachromatic(self):
        """
        This hit fails to identify that the extra chloride comes from x1382.
        The human chose

        * the methylpyridine off x0107
        * the benzene-pyridine off x0434, but wanted an amide not a ureido
        * the link between the rings as amide x0678
        * x0995 is a red herring
        * benzene with a chroride from x1382


        :return:
        """
        MProVictor.quick_renanimation = True
        # ,'x2646'
        Victor.monster_throw_on_discard = True
        victor = MProVictor.from_hit_codes(smiles='Cc1ccncc1NC(=O)Cc1cccc(Cl)c1',
                                           # hit_codes=['x0107','x0434','x0678','x0995','x1382'],
                                           hit_codes=['x0107', 'x0434', 'x1382'],
                                           long_name='TRY-UNI-714a760b-6')
        self.assertEqual(victor.error, '', victor.error)
        self.assertIsNotNone(victor.minimised_mol, 'Failed minimisation')
        actual = MProVictor.get_mol('x2646')
        victor.make_pse(extra_mols=[actual])
        rmsd = victor.validate(reference_mol=actual)
        self.assertLess(rmsd, 1, f'The RMSD is large...')
        # self.assertIn('x1382', victor.monster.matched)
        # self.assertIn('x0995', victor.monster.unmatched) # red herring


# ======================================================================================================================


class RectifierTester(unittest.TestCase):
    def test_rectifier(self):
        # name: [before, after]
        chemdex = {'phenylnaphthalene': ('c1ccc2ccccc2c1(c3ccccc3)', 'c1ccc(-c2cccc3ccccc23)cc1'),
                   'benzo-azetine': ('C12CCCCC1CC2', 'C1CCC2CCCC2C1'),
                   # 'conjoined': ('C1C2CCC2C1', 'C1CCCCC1'), # bridged hexane
                   'allene': ('C=C=C', 'C=CC'),
                   'benzo-cyclopronane': ('C12CCCCC1C2', 'C1CCC2CCCC2C1'),
                   # 'norbornane': ('C1CC2CCC1C2', 'C1CC2CCC1C2'),
                   'mixed ring': ('c1cccc2c1CCCC2', 'c1ccc2c(c1)CCCC2'),
                   'mixed ring': ('C1CCCc2c1cccc2', 'c1ccc2c(c1)CCCC2'),
                   }

        for name in chemdex:
            before, after = chemdex[name]
            mol = Chem.MolFromSmiles(before)
            mol.SetProp('_Name', name)
            AllChem.EmbedMolecule(mol)
            recto = Rectifier(mol).fix()
            gotten = Chem.MolToSmiles(recto.mol)
            self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after}) from {before}')

    def test_cyclopentine(self):
        # aromatic cyclopent-ine -> cyclopentadiene
        name = 'cyclopentine'
        mol = Chem.MolFromSmiles('[nH]1cccc1')
        mol.SetProp('_Name', name)
        AllChem.EmbedMolecule(mol)
        mod = Chem.RWMol(mol)
        mod.GetAtomWithIdx(0).SetAtomicNum(6)
        mol = mod.GetMol()
        recto = Rectifier(mol, atoms_in_bridge_cutoff=3).fix()
        gotten = Chem.MolToSmiles(recto.mol)
        after = 'C1=CCC=C1'
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_bad_ring(self):
        name = 'bad ring'
        after = 'c1ccc2c(c1)CCCC2'
        mol = Chem.MolFromSmiles(after)
        mol.SetProp('_Name', name)
        mol.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.SINGLE)
        before = Chem.MolToSmiles(mol)
        recto = Rectifier(mol).fix()
        gotten = Chem.MolToSmiles(recto.mol)
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_bad_ring2(self):
        name = 'bad ring2'
        before = 'c1ccc2c(c1)CCCC2'
        after = 'c1ccc2ccccc2c1'
        mol = Chem.MolFromSmiles(before)
        mol.SetProp('_Name', name)
        mol.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.SINGLE)
        mol.GetBondBetweenAtoms(6, 7).SetBondType(Chem.BondType.AROMATIC)
        before = Chem.MolToSmiles(mol)
        recto = Rectifier(mol).fix()
        gotten = Chem.MolToSmiles(recto.mol)
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')


class RingTestsVictor(unittest.TestCase):

    def test_orthomethyltoluene(self):
        name = 'orthomethyltoluene'
        after = 'Cc1cccc(C)c1'
        template = os.path.join(MProVictor.get_mpro_path(), 'template.pdb')
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        rototoluene = Chem.MolFromMolFile('test_mols/rototoluene.mol')
        rototoluene.SetProp('_Name', 'rototoluene')
        victor = Victor.combine(hits=[toluene, rototoluene],
                           pdb_filename=template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')
        self.assertEqual(victor.error, '', victor.error)
        gotten = Chem.MolToSmiles(Chem.RemoveHs(victor.minimised_mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_peridimethylnaphthalene(self):
        name = 'peridimethylnaphthalene'
        after = 'Cc1cccc2cccc(C)c12'
        template = os.path.join(MProVictor.get_mpro_path(), 'template.pdb')
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneF = Chem.MolFromMolFile('test_mols/transtoluene.mol')
        transtolueneF.SetProp('_Name', 'transtoluene-fuse')
        victor = Victor.combine(hits=[toluene, transtolueneF],
                           pdb_filename=template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')
        self.assertEqual(victor.error, '', victor.error)
        gotten = Chem.MolToSmiles(Chem.RemoveHs(victor.minimised_mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_spirodituluene(self):
        name = 'spirodituluene'
        after = ('C[C@@H]1C=CC[C@]2(C=C[C@H](C)C=C2)C1', 'C[C@@H]1C=CC[C@]2(C=C[C@@H](C)C=C2)C1')
        template = os.path.join(MProVictor.get_mpro_path(), 'template.pdb')
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneS = Chem.MolFromMolFile('test_mols/transtoluene2.mol')
        transtolueneS.SetProp('_Name', 'transtoluene-spiro')
        # cmd.rotate('z', -90, 'rototoluene', camera=0)
        victor = Victor.combine(hits=[toluene, transtolueneS],
                           pdb_filename=template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')
        self.assertEqual(victor.error, '', victor.error)
        gotten = Chem.MolToSmiles(Chem.RemoveHs(victor.minimised_mol))
        self.assertIn(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_phenylene(self):
        # make carboxy and amide benzenes that overlap so that the end result is a phenylene where one ring is oxazine
        before = 'c3c1cccc2\C(=O)O/C(-N)c(c12)cc3'
        conjoined = Chem.MolFromSmiles(before)
        AllChem.EmbedMolecule(conjoined)
        bonds = [conjoined.GetBondBetweenAtoms(0, 1).GetIdx(),
                 conjoined.GetBondBetweenAtoms(12, 11).GetIdx(),
                 conjoined.GetBondBetweenAtoms(8, 9).GetIdx()]
        fragged = Chem.FragmentOnBonds(conjoined, bonds, addDummies=False)
        fore = Chem.GetMolFrags(fragged, asMols=True, sanitizeFrags=False)[1]
        Chem.SanitizeMol(fore)
        bonds = [conjoined.GetBondBetweenAtoms(2, 1).GetIdx(),
                 conjoined.GetBondBetweenAtoms(12, 5).GetIdx(),
                 conjoined.GetBondBetweenAtoms(8, 6).GetIdx()]
        fragged = Chem.FragmentOnBonds(conjoined, bonds, addDummies=False)
        aft = Chem.GetMolFrags(fragged, asMols=True, sanitizeFrags=False)[0]
        Chem.SanitizeMol(aft)
        # merge them
        mol = Monster([fore, aft]).merge().positioned_mol
        after = Chem.MolToSmiles(mol)
        self.assertEqual(before, after)


# ----------------------------------------------------------------------------------------------------------------------


class Internals(unittest.TestCase):
    def test_triangle(self):
        """
        Test triangle prevention.
        """
        probanda = Chem.MolFromSmiles('CCC')
        AllChem.EmbedMolecule(probanda)
        monster = Monster([probanda])
        zeroth = probanda.GetAtomWithIdx(0)
        second = probanda.GetAtomWithIdx(2)
        # 0 - 2 would make a triangle
        self.assertTrue(
            monster._is_would_be_triangle(zeroth, second))  # connecting zeroth and second would make a triangle
        self.assertEqual(monster._get_triangle(zeroth, second), 1)  # connecting 0+2, would make 1 the vertex.

    def make_mol(self, smiles: str) -> Chem.Mol:
        mol = Chem.MolFromSmiles(smiles)
        dummies = mol.GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0))
        for dummy in dummies:
            dummy.SetAtomicNum(6)  # carbon is 6 not 12!
        AllChem.EmbedMolecule(mol)
        for dummy in dummies:
            dummy.SetAtomicNum(0)
        return mol

    def make_pair_by_split(self, conjoined: Chem.Mol, atom_idx: int) -> Tuple[Chem.Mol]:
        # make overlapping mols by getting a single molecule, and split it
        # this gives more control over Chem.rdMolAlign.AlignMol as this may overlap other atoms.
        # negative weights does not work...
        # fore
        bond = conjoined.GetBondBetweenAtoms(atom_idx, atom_idx + 1)
        fragged = Chem.FragmentOnBonds(conjoined, [bond.GetIdx()], addDummies=False)
        fore = Chem.GetMolFrags(fragged, asMols=True)[0]
        bond = conjoined.GetBondBetweenAtoms(atom_idx - 1, atom_idx)
        fragged = Chem.FragmentOnBonds(conjoined, [bond.GetIdx()], addDummies=False)
        aft = Chem.GetMolFrags(fragged, asMols=True)[1]
        return fore, aft

    def test_merge_on_same_dummy(self):
        conjoined = self.make_mol('O=C(O)CSCC#N')
        acetyl, nitrile = self.make_pair_by_split(conjoined, 4)
        # merge
        monster = Monster([acetyl, nitrile])
        merger = monster.merge_hits()
        dummies = merger.GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0))
        for dummy in dummies:
            self.assertEqual(len(dummy.GetNeighbors()), 1)
        self.assertEqual(len(Chem.GetMolFrags(merger)), 1)


# ----------------------------------------------------------------------------------------------------------------------
class UnresolvedProblems(unittest.TestCase):
    def test_recto_fail_A(self):
        """Not too sure why this fails. I think it is the alphatic - ring merger"""
        MProVictor.monster_throw_on_discard = True
        victor = MProVictor.combine_codes(hit_codes=['x11612', 'x11475'])
        self.assertEqual(victor.error, '', victor.error)

    def test_supplementary1_to_recto_fail_A(self):
        """
        This was ment to test the above, but it works fine.
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
        monster = Monster(hits=[toluene, chlorobutane]).merge(keep_all=False)
        # ======
        self.assertEqual(Chem.MolToSmiles(monster.positioned_mol), Chem.MolToSmiles(chlorotoluene))  # CC(Cl)CCc1ccccc1

    def test_supplementary2_to_recto_fail_A(self):
        """
        This was meant to test as above. It also works fine.
        :return:
        """
        #
        methylchlorotoluene = Chem.MolFromSmiles('Cc1c(Cl)c(C)ccc1')
        AllChem.EmbedMolecule(methylchlorotoluene)
        methylchlorotoluene.SetProp('_Name', 'methylchlorotoluene')
        #
        methyltoluene = Chem.RWMol(methylchlorotoluene)
        methyltoluene.RemoveAtom(3)
        Chem.SanitizeMol(methyltoluene)
        methyltoluene.SetProp('_Name', 'methyltoluene')
        #
        chloropentane = Chem.RWMol(methylchlorotoluene)
        for n in range(chloropentane.GetNumAtoms() - 1, 5, -1):
            chloropentane.RemoveAtom(n)
        for atom in chloropentane.GetAtoms():
            atom.SetIsAromatic(False)
        for bond in chloropentane.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
        Chem.SanitizeMol(chloropentane)
        chloropentane.SetProp('_Name', '2-chloropentane')
        #
        monster = Monster(hits=[methyltoluene, chloropentane]).merge(keep_all=False)
        # ======
        self.assertEqual(Chem.MolToSmiles(monster.positioned_mol), Chem.MolToSmiles(methylchlorotoluene))

# Todo: add a class to test missing modules.
