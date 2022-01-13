import logging
import tempfile
import unittest, os
# ======================================================================================================================
import pyrosetta

pyrosetta.init(
    extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')
# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem

from fragmenstein import Monster, Victor, Igor
from fragmenstein.mpro import MProVictor
from typing import *
import numpy as np


# ======================================================================================================================


class MProPlaceTester(unittest.TestCase):

    def test_easy(self):
        """
        To a **human** this looks easy. x0692 is a red herring and the other two make the molecule.
        As currently written this will fail.

        :return:
        """
        # PAU-UNI-52c0427f-1
        MProVictor.quick_reanimation = True
        victor = MProVictor.from_hit_codes(hit_codes=['x0692', 'x0305', 'x1249'])
        victor.place(smiles='CCNc1ncc(C#N)cc1CN1CCN(C(=O)C*)CC1', long_name='2_ACL')
        self.assertEqual(victor.error_msg, '', victor.error_msg)
        self.assertIsNotNone(victor.minimized_mol, 'Failed minimisation')
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
        MProVictor.quick_reanimation = True
        victor = MProVictor.from_hit_codes(hit_codes=['x0434', 'x0540'])
        victor.place(smiles='*CCC(=O)N1CC(CCN(C(=O)Nc2c(C)ncc(C)c2CCN2CCOCC2)c2cc(C)ccn2)C1',
                     long_name='AGN-NEW-5f0-1_ACR1')
        self.assertEqual(str(victor.error_msg), '', str(victor.error_msg))
        self.assertIsNotNone(victor.minimized_mol, 'Failed minimisation')
        self.assertEqual(len(victor.monster.unmatched), 0,
                         f'Both were correct but {victor.monster.unmatched} was discarded')
        victor.make_pse()

    def test_incorrect(self):
        """
        This case has two hits that are not ideal. One, x0995, totally barney.
        These will be rejected.

        :return:
        """
        MProVictor.quick_reanimation = True
        victor = MProVictor.from_hit_codes(hit_codes='x0692,x0770,x0995'.split(','))
        victor.place(smiles='*C(=N)CN1CCN(Cc2ccc(-c3cc(CC)ncn3)c(F)c2)CC1',
                     long_name='BEN-VAN-c98-4')
        self.assertEqual(victor.error_msg, '', victor.error_msg)
        self.assertIsNotNone(victor.minimized_mol, 'Failed minimisation')
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
        MProVictor.quick_reanimation = True
        # ,'x2646'
        Victor.monster_throw_on_discard = True
        victor = MProVictor.from_hit_codes(# hit_codes=['x0107','x0434','x0678','x0995','x1382'],
                                           hit_codes=['x0107', 'x0434', 'x1382'])
        victor.place(smiles='Cc1ccncc1NC(=O)Cc1cccc(Cl)c1',
                     long_name='TRY-UNI-714a760b-6')
        self.assertEqual(victor.error_msg, '', victor.error_msg)
        self.assertIsNotNone(victor.minimized_mol, 'Failed minimisation')
        actual = MProVictor.get_mol('x2646')
        victor.make_pse(extra_mols=[actual])
        rmsd = victor.validate(reference_mol=actual)
        self.assertLess(rmsd, 1, f'The RMSD is large...')
        # self.assertIn('x1382', victor.monster.matched)
        # self.assertIn('x0995', victor.monster.unmatched) # red herring


# ======================================================================================================================

class VictorCombineTests(unittest.TestCase):
    def test_noncovalent(self):
        MProVictor.quick_reanimation = False
        victor = MProVictor.from_hit_codes(hit_codes=['x0305', 'x1249'])
        victor.combine()
        self.assertLess(victor.mrmsd.mrmsd, 1, f'RMSD great that one ({victor.mrmsd.mrmsd})')
        self.assertLess(victor.ddG, -1, f'ddG {victor.ddG}')

    def test_covalent(self):
        MProVictor.quick_reanimation = False
        victor = MProVictor.from_hit_codes(hit_codes=['x0692', 'x0305', 'x1249'])
        victor.combine()
        self.assertLess(victor.mrmsd.mrmsd, 1, f'RMSD great that one ({victor.mrmsd.mrmsd})')
        self.assertLess(victor.ddG, -1, f'ddG {victor.ddG}')

class MonsterCombineTests(unittest.TestCase):
    def test_phenylene(self):
        # make carboxy and amide benzenes that overlap so that the end result is a phenylene where one ring is oxazine
        conjoined = Chem.MolFromSmiles('c3c1cccc2\C(=O)O/C(-N)c(c12)cc3')
        before = Chem.MolToSmiles(conjoined) # structure from wiki is not canonical
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
        mol = Monster([fore, aft]).combine().positioned_mol
        after = Chem.MolToSmiles(mol)
        self.assertEqual(before, after)

    def test_orthomethyltoluene(self):
        name = 'orthomethyltoluene'
        after = 'Cc1cccc(C)c1'
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        rototoluene = Chem.MolFromMolFile('test_mols/rototoluene.mol')
        rototoluene.SetProp('_Name', 'rototoluene')
        mol = Monster(hits=[toluene, rototoluene]).combine(keep_all=True).positioned_mol
        gotten = Chem.MolToSmiles(Chem.RemoveHs(mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_peridimethylnaphthalene(self):
        name = 'peridimethylnaphthalene'
        after = 'Cc1cccc2cccc(C)c12'
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneF = Chem.MolFromMolFile('test_mols/transtoluene.mol')
        transtolueneF.SetProp('_Name', 'transtoluene-fuse')
        mol=Monster(hits=[toluene, transtolueneF]).combine(keep_all=True).positioned_mol
        gotten = Chem.MolToSmiles(Chem.RemoveHs(mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_spirodituluene(self):
        name = 'spirodituluene'
        after = ('CC1C=CC2(C=C1)CC=CC(C)C2',
                 'C[C@@H]1C=CC[C@]2(C=C[C@H](C)C=C2)C1',
                 'C[C@@H]1C=CC[C@]2(C=C[C@@H](C)C=C2)C1',
                 'C[C@H]1C=CC[C@]2(C=C[C@H](C)C=C2)C1',
                 'C[C@H]1C=CC[C@]2(C=C[C@@H](C)C=C2)C1')
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneS = Chem.MolFromMolFile('test_mols/transtoluene2.mol')
        transtolueneS.SetProp('_Name', 'transtoluene-spiro')
        # cmd.rotate('z', -90, 'rototoluene', camera=0)
        mol=Monster(hits=[toluene, transtolueneS]).combine(keep_all=True).positioned_mol
        gotten = Chem.MolToSmiles(Chem.RemoveHs(mol))
        self.assertIn(gotten, after, f'{name} failed {gotten} (expected {after})')


class MonsterPlaceTests(unittest.TestCase):

    def test_sample_new_conformation(self):
        smiles = "C1C2C(C=C(C=2)C(C2C=CC=C2)CNOC)C=CC=1"
        mol = Chem.MolFromSmiles(smiles)
        AllChem.EmbedMolecule(mol, randomSeed=131)
        frag_mol = Chem.FragmentOnBonds(mol, [5, 11], addDummies=False)
        hits = Chem.GetMolFrags(frag_mol, asMols=True)[:2]
        # from matplotlib import pyplot as plt
        # from rdkit.Chem import Draw
        # plt.imshow(Draw.MolsToGridImage(hits)); plt.show()
        ori_monster =  Monster(hits=hits, random_seed=131)
        ori_monster.place_smiles(smiles)
        # sample_new_conformation
        seeds = [121, 23421, 1]
        for se1 in seeds:
            mol = ori_monster.sample_new_conformation(random_seed=se1)
            coords1 = mol.GetConformer().GetPositions()
            for se2 in seeds:
                mol = ori_monster.sample_new_conformation(random_seed=se2)
                coords2 = mol.GetConformer().GetPositions()
                if se1==se2:
                    self.assertAlmostEqual( np.sum(np.abs(coords1-coords2)), 0)
                else:
                    self.assertTrue( np.sum(np.abs(coords1-coords2)) > 3 )

    def test_random_seed(self):
        smiles = "C1C2C(C=C(C=2)C(C2C=CC=C2)CNOC)C=CC=1"
        mol = Chem.MolFromSmiles(smiles)
        AllChem.EmbedMolecule(mol, randomSeed=131)
        frag_mol = Chem.FragmentOnBonds(mol, [5, 11], addDummies=False)
        hits = Chem.GetMolFrags(frag_mol, asMols=True)[:2]
        # from matplotlib import pyplot as plt
        # from rdkit.Chem import Draw
        # plt.imshow(Draw.MolsToGridImage(hits)); plt.show()
        seeds = [121, 23421, 1]
        for se1 in seeds:
            mol = Monster(hits=hits, random_seed=se1).place_smiles(smiles).positioned_mol
            coords1 = mol.GetConformer().GetPositions()
            for se2 in seeds:
                mol = Monster(hits=hits, random_seed=se2).place_smiles(smiles).positioned_mol
                coords2 = mol.GetConformer().GetPositions()
                if se1==se2:
                    self.assertAlmostEqual( np.sum(np.abs(coords1-coords2)), 0)
                else:
                    self.assertTrue( np.sum(np.abs(coords1-coords2)) > 3 )

class VictorPlaceTests(unittest.TestCase):
    def test_random_seed(self):

        to_place = Chem.MolFromMolFile('test_mols/placed_example1.mol')
        pdb_filename = 'test_mols/apo_example1.pdb'
        smiles = Chem.MolToSmiles(to_place)
        hits = [ Chem.MolFromMolFile(os.path.join('test_mols', basename)) for basename in ["x0032_0A.mol"]] #, "x0103_0A.mol"]]
        seeds = [121, 23421, 1]
        # Victor.enable_stdout(level=logging.ERROR)
        Victor.monster_mmff_minisation = False
        import random

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir= "/home/sanchezg/tmp/pruebaFragmenstein"
            for se1 in seeds:
                from pyrosetta.rosetta.basic.random import determine_random_number_seed, init_random_generators
                init_random_generators(se1, "mt19937")
                random.seed(se1)
                np.random.seed(se1)
                Victor.work_path = os.path.join(tmpdir, "out_1_%d" % (se1))
                mol = Victor(hits=hits, pdb_filename=pdb_filename, random_seed=se1).place(smiles).minimized_mol
                coords1 = mol.GetConformer().GetPositions()
                for se2 in seeds:
                    Victor.work_path = os.path.join(tmpdir, "out_2_%d" % (se2))
                    # pyrosetta.init(extra_options=f'-constant_seed true -jran 987')
                    init_random_generators(se2, "mt19937")
                    random.seed(se2)
                    np.random.seed(se2)
                    mol = Victor(hits=hits, pdb_filename=pdb_filename, random_seed=se2).place(smiles).minimized_mol
                    coords2 = mol.GetConformer().GetPositions()
                    if se1==se2:
                        self.assertAlmostEqual( np.sum(np.abs(coords1-coords2)), 0)
                    else:
                        print(se1, se2, np.sum(np.abs(coords1-coords2)))
                        self.assertTrue( np.sum(np.abs(coords1-coords2)) > 3 )

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

    def make_pair_by_split(self, conjoined: Chem.Mol, atom_idx: int) -> Tuple[Chem.Mol,Chem.Mol]:
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
        merger = monster.simply_merge_hits()
        dummies = merger.GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0))
        for dummy in dummies:
            self.assertEqual(len(dummy.GetNeighbors()), 1)
        self.assertEqual(len(Chem.GetMolFrags(merger)), 1)

    def translate(self, mol, x=0, y=0, z=0):
        """
        Translates the molecule in place

        :param mol:
        :param x:
        :param y:
        :param z:
        :return:
        """
        translation = np.array([[1, 0, 0, x],
                                [0, 1, 0, y],
                                [0, 0, 1, z],
                                [0, 0, 0, 1]], dtype=np.double)
        AllChem.TransformConformer(mol.GetConformer(0), translation)

    def test_join(self):
        wanted = 'c1ccc(Oc2ccccc2)cc1'
        benzene = Chem.MolFromSmiles('c1ccccc1')
        AllChem.EmbedMolecule(benzene)
        moved = Chem.Mol(benzene)
        self.translate(moved, x=5)
        found = Chem.MolToSmiles(Monster([benzene, moved]).combine().positioned_mol)
        self.assertEqual(wanted, found, 'The joining differs')

    def test_distance(self):
        methane = Chem.MolFromSmiles('C')
        AllChem.EmbedMolecule(methane)
        ammonia = Chem.MolFromSmiles('N')
        AllChem.EmbedMolecule(ammonia)
        self.translate(ammonia, x=3)
        monster = Monster([methane, ammonia])
        monster.combine(keep_all=False, joining_cutoff=2)
        self.assertEqual(1, len(monster.unmatched), 'discard error')
        monster.combine(keep_all=True, joining_cutoff=5)
        self.assertEqual(0, len(monster.unmatched), 'discard error')
        try:
            monster.combine(keep_all=True, joining_cutoff=2)
            self.fail('should have raised a connection error')
        except ConnectionError as error:
            pass
# ----------------------------------------------------------------------------------------------------------------------
class UnresolvedProblems(unittest.TestCase):
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
        self.assertEqual(Chem.MolToSmiles(chlorotoluene), Chem.MolToSmiles(monster.positioned_mol))  # CC(Cl)CCc1ccccc1

    def test_supplementary2_to_recto_fail_A(self):
        """
        This was meant to test as above.
        It mergers xylene with chloropentane
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
        self.assertEqual(Chem.MolToSmiles(chloroxylene), Chem.MolToSmiles(monster.positioned_mol))

# Todo: add a class to test missing modules.
