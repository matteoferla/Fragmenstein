import logging
import tempfile
import unittest, os
# ======================================================================================================================
from multiprocessing import Process

import pyrosetta

# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem

from fragmenstein import Monster, Victor, Igor, mpro_data, Walton
from fragmenstein.mpro import MProVictor
from typing import *
import numpy as np

# ======================================================================================================================

import os
test_mols_folder = os.path.join(os.path.dirname(__file__), 'test_mols')

class MProPlaceTester(unittest.TestCase):

    def setUp(self):
        Igor.init_pyrosetta()

    def untest_red_herring(self):  # without the test_ word this will not run.
        """
        To a **human** this looks easy. x0692 is a red herring and the other two make the molecule.
        As currently written this will fail.

        See red herring test notes in documentation.

        :return:
        """
        # PAU-UNI-52c0427f-1
        MProVictor.quick_reanimation = True
        victor = MProVictor.from_hit_codes(hit_codes=['x0692', 'x0305', 'x1249'])
        victor.place(smiles='CCNc1ncc(C#N)cc1CN1CCN(C(=O)C*)CC1', long_name='2_ACL')
        self.assertEqual(victor.error_msg, '', victor.error_msg)
        self.assertIsNotNone(victor.minimized_mol, 'Failed minimisation')
        msg = f'x1249 is the red herring, prediction: {victor.monster.unmatched} discarded, ' + \
              f'while x0305 and x0692 are the true inspirations. kept: {victor.monster.matched}'
        self.assertIn('x1249', victor.monster.unmatched, msg)
        self.assertIn('x0305', victor.monster.matched, msg)

        victor.make_pse()

    def untest_nasty(self):  # without the test_ word this will not run.
        """
        The human suggested a lot of novel groups.
        'x0540' is a really odd inspiration. Three atoms are conserved. the rest aren't.

        Like the red herring this is impossible.
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
        The human chose these.

        * the methylpyridine off x0107
        * the benzene-pyridine off x0434, but wanted an amide not a ureido
        * the link between the rings as amide x0678
        * x0995 is a red herring
        * benzene with a chroride from x1382
        """
        MProVictor.quick_reanimation = True
        # ,'x2646'
        Victor.monster_throw_on_discard = True
        victor = MProVictor.from_hit_codes(  # hit_codes=['x0107','x0434','x0678','x0995','x1382'],
            hit_codes=['x0107', 'x0434', 'x1382'])
        victor.place(smiles='Cc1ccncc1NC(=O)Cc1cccc(Cl)c1',
                     long_name='TRY-UNI-714a760b-6')
        self.assertEqual(victor.error_msg, '', victor.error_msg)
        self.assertIsNotNone(victor.minimized_mol, 'Failed minimisation')
        actual = mpro_data.get_mol('x2646')
        victor.make_pse(extra_mols=[actual])
        validation: Dict[str, float] = victor.validate(reference_mol=actual)
        rmsd = validation['reference2minimized_rmsd']
        self.assertLess(rmsd, 2, f'The RMSD is large...')
        # self.assertIn('x1382', victor.monster.matched)
        # self.assertIn('x0995', victor.monster.unmatched) # red herring


# ======================================================================================================================

class VictorCombineTests(unittest.TestCase):

    def setUp(self):
        Igor.init_pyrosetta()

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
        self.assertLess(victor.mrmsd.mrmsd, 1.2, f'RMSD great that one ({victor.mrmsd.mrmsd})')
        self.assertLess(victor.ddG, -1, f'ddG {victor.ddG}')


class MonsterCombineTests(unittest.TestCase):
    def test_phenylene(self):
        # make carboxy and amide benzenes that overlap so that the end result is a phenylene where one ring is oxazine
        conjoined = Chem.MolFromSmiles('c3c1cccc2\C(=O)O/C(-N)c(c12)cc3')
        conjoined = next(AllChem.EnumerateStereoisomers(conjoined))
        before = Chem.MolToSmiles(conjoined)  # structure from wiki is not canonical
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
        after = Chem.MolToSmiles(Chem.RemoveHs(mol))
        self.assertEqual(before, after)

    def test_orthomethyltoluene(self):
        name = 'orthomethyltoluene'
        after = 'Cc1cccc(C)c1'
        toluene = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'toluene.mol'))
        toluene.SetProp('_Name', 'toluene')
        rototoluene = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'rototoluene.mol'))
        rototoluene.SetProp('_Name', 'rototoluene')
        mol = Monster(hits=[toluene, rototoluene]).combine(keep_all=True).positioned_mol
        gotten = Chem.MolToSmiles(Chem.RemoveHs(mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_peridimethylnaphthalene(self):
        name = 'peridimethylnaphthalene'
        after = 'Cc1cccc2cccc(C)c12'
        toluene = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'toluene.mol'))
        toluene.SetProp('_Name', 'toluene')
        transtolueneF = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'transtoluene.mol'))
        transtolueneF.SetProp('_Name', 'transtoluene-fuse')
        mol = Monster(hits=[toluene, transtolueneF]).combine(keep_all=True).positioned_mol
        gotten = Chem.MolToSmiles(Chem.RemoveHs(mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_spirodituluene(self):
        name = 'spirodituluene'
        after = ('CC1CCC2(CCCC(C)C2)CC1',
                 'CC1C=CC2(C=C1)CC=CC(C)C2',
                 'C[C@@H]1C=CC[C@]2(C=C[C@H](C)C=C2)C1',
                 'C[C@@H]1C=CC[C@]2(C=C[C@@H](C)C=C2)C1',
                 'C[C@H]1C=CC[C@]2(C=C[C@H](C)C=C2)C1',
                 'C[C@H]1C=CC[C@]2(C=C[C@@H](C)C=C2)C1')
        toluene = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'toluene.mol'))
        toluene.SetProp('_Name', 'toluene')
        transtolueneS = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'transtoluene2.mol'))
        transtolueneS.SetProp('_Name', 'transtoluene-spiro')
        # cmd.rotate('z', -90, 'rototoluene', camera=0)  # this predates capt. Robert Walton
        mol = Monster(hits=[toluene, transtolueneS]).combine(keep_all=True).positioned_mol
        gotten = Chem.MolToSmiles(Chem.RemoveAllHs(mol))
        self.assertIn(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_real_merger(self):
        # Victor.enable_stdout(logging.DEBUG)
        x0138 = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'mac-x0138.mol'))
        x0398 = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'mac-x0398.mol'))
        monster = Monster([x0398, x0138])
        monster.combine()
        self.assertEqual('Nc1nc2c3c(c(O)cc(N)c3n1)C(O)=N2',
                         Chem.MolToSmiles(Chem.RemoveHs(monster.positioned_mol)),
                         )


class MonsterPlaceTests(unittest.TestCase):

    def get_5SB7_mols(self) -> List[Chem.Mol]:
        """
        In https://onlinelibrary.wiley.com/doi/10.1002/anie.202204052
        There is a lovely test case.
        Namely the followup bold4 is a merger of F36 and F04 hits.
        However... F04 (hello Hantzsch thiazole synthesis product) presents
        the problematic case that the thiazole can be flipped either way
        upon investigation of the crystal map
        and the followup bold4 is indeed a flipped merger.
        Therefore, without a custom_map the followup is incorrect.
        """
        with Chem.SDMolSupplier('test_mols/5SB7_mergers.sdf') as reader:
            return list(reader)


    def test_sample_new_conformation(self):
        smiles = "C1C2C(C=C(C=2)C(C2C=CC=C2)CNOC)C=CC=1"
        mol = Chem.MolFromSmiles(smiles)
        AllChem.EmbedMolecule(mol, randomSeed=131)
        frag_mol = Chem.FragmentOnBonds(mol, [5, 11], addDummies=False)
        hits = Chem.GetMolFrags(frag_mol, asMols=True)[:2]
        # from matplotlib import pyplot as plt
        # from rdkit.Chem import Draw
        # plt.imshow(Draw.MolsToGridImage(hits)); plt.show()
        ori_monster = Monster(hits=hits, random_seed=131)
        ori_monster.place_smiles(smiles)
        # sample_new_conformation
        seeds = [121, 23421, 1]
        for se1 in seeds:
            mol = ori_monster.sample_new_conformation(random_seed=se1)
            coords1 = mol.GetConformer().GetPositions()
            for se2 in seeds:
                mol = ori_monster.sample_new_conformation(random_seed=se2)
                coords2 = mol.GetConformer().GetPositions()
                if se1 == se2:
                    self.assertAlmostEqual(np.sum(np.abs(coords1 - coords2)), 0)
                else:
                    self.assertTrue(np.sum(np.abs(coords1 - coords2)) > 3)

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
                if se1 == se2:
                    self.assertAlmostEqual(np.sum(np.abs(coords1 - coords2)), 0)
                else:
                    self.assertTrue(np.sum(np.abs(coords1 - coords2)) > 3)

    def test_flipped_lactam(self):
        """
        Given a benzo + 7-membered lactam map a mol with the amide flipped
        """
        hit_F584 = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'F584.mol'))

        monster = Monster(hits=[hit_F584, ])
        monster.place_smiles(smiles='COc1cccc2C(=O)NCCCc12', long_name='flipped_F584')
        # monster.show_comparison()
        self.assertEqual(len(monster.get_mcs_mappings(followup=monster.initial_mol, hit=monster.hits[0])[0][0]), 13)

    def test_forced_flipped_lactam(self):

        """
        Given a benzo + 7-membered lactam map a mol with the amide flipped as the previous test
        force the ketones to match
        """
        hit_F584 = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'F584.mol'))
        monster = Monster(hits=[hit_F584, ])
        monster.place_smiles(smiles='COc1cccc2C(=O)NCCCc12',
                             long_name='flipped_F584',
                             custom_map={'hit_F584': {13: 8}},
                             # merging_mode='expansion', # works on both expansion and no blend modes
                             )
        # monster.draw_nicely(monster.hits[0])
        # monster.draw_nicely(monster.initial_mol)
        # monster.show_comparison()
        # the [0][0][13] is because of:
        # Tuple[List[Dict[int, int]], ExtendedFMCSMode]
        self.assertIn(13, monster.get_mcs_mappings(followup=monster.initial_mol, hit=monster.hits[0])[0][0])
        self.assertEqual(monster.get_mcs_mappings(followup=monster.initial_mol, hit=monster.hits[0])[0][0][13], 8)
        #self.assertEqual(len(monster.get_mcs_mappings(followup=monster.initial_mol, hit=monster.hits[0])[0][0]), 13)

    def test_thiazole_flip(self):
        """
        See ``.get_5SB7_mols`` for details.

        Potentially the number of specified atoms is excessive,
        but as a human I just wrote down the atoms in the ring and was done with it.
        """
        mols = self.get_5SB7_mols()
        monster = Monster([mols[0]])
        monster.place(Chem.Mol(mols[0]),
                      merging_mode='expansion',
                      # custom_map={'F36': {1:7}, 'F04': {4:7}}
                      custom_map={'F04': {-1: 4,  # no amine
                                          4: -2,  # no amine
                                          12: 12,
                                          6: 13,
                                          13: 6,
                                          15: 14,
                                          14: 15}}
                      )
        self.assertEqual(len(list(filter(len, monster.origin_from_mol(monster.positioned_mol)))), 15)
        # the amine is banned in the map and the ring is flipped.

    def test_thiazole_followup(self):
        """
        The followup compound is from 5SB3
        and is 'bold4' in mols.
        """
        mols = self.get_5SB7_mols()
        monster = Monster(mols[:2])
        monster.place(Chem.Mol(mols[3]),
                      merging_mode='expansion',
                      custom_map={'F36': {1: 7},
                                  'F04': {4: -1,  # no amine
                                          12: 13,  # root to Ph
                                          13: 5,
                                          6: 14,
                                          }
                                  }
                      )
        # monster.show_comparison()
        # monster.to_nglview(True)
        # the amine at 7 is banned in F04.
        self.assertEqual(monster.origin_from_mol()[7], ['F36.1'])
        self.assertEqual(len(list(filter(len, monster.origin_from_mol(monster.positioned_mol)))), 22)


class MultivictorPlaceTests(unittest.TestCase):

    def setUp(self):
        Igor.init_pyrosetta()

    def test_multivictor(self):
        from fragmenstein import MultiVictorPlacement
        to_place = Chem.MolFromMolFile(os.path.join(test_mols_folder, 'placed_example1.mol'))
        pdb_filename = os.path.join(test_mols_folder, 'apo_example1.pdb')
        smiles = Chem.MolToSmiles(to_place)
        hits = [Chem.MolFromMolFile(os.path.join(test_mols_folder, basename)) for basename in
                ["x0032_0A.mol"]]  # , "x0103_0A.mol"]]
        # Victor.enable_stdout(level=logging.ERROR)
        with tempfile.TemporaryDirectory() as tmpdir:
            Victor.work_path = os.path.join(tmpdir, "multivictor_out")
            mv = MultiVictorPlacement(hits=hits, pdb_filename=pdb_filename)
            mv.place(smiles, number_runs=4)
            # print(mv.retrieve_best_victor())
            # print(mv.retrieve_scores())
            self.assertLess(mv.retrieve_scores()[0], -7)

class WaltonTests(unittest.TestCase):

    def smiles_assertEqual(self, a, b):
        """
        helper method to test equality of smiles
        """
        if isinstance(a, str):
            expected_smiles = a
            mol = b
        else:
            expected_smiles = b
            mol = a
        obtained_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
        self.assertEqual(expected_smiles, obtained_smiles)

    def test_align_map(self):
        demo = Walton.from_smiles(resorcinol='c1ccc(O)cc1O', eugenol='Oc1ccc(cc1OC)CC=C')  # create instance
        demo.align_by_map({(0, 1): {4: 0, 3: 1, 2: 2}})  # align molecules by atom indices
        demo()  # merge (Fragmenstein's Monster)
        self.smiles_assertEqual('C=CCc1ccc(O)c(OC)c1O', demo.merged)

    def test_polygon(self):
        benzene = Walton.create_polygon(n=6, bond_type=Chem.BondType.AROMATIC)
        self.smiles_assertEqual('c1ccccc1', benzene)

    def test_ring_on_plane(self):
        demo = Walton.from_smiles(furan='c1ccco1')
        demo.ring_on_plane(mol_idx=0)  # flat on xy
        demo.atom_on_axis(mol_idx=0, atom_idx=4, axis='x')  # ox on axis
        # the atom_idx can be actually be passed a Geometry.Point3D:
        demo.atom_on_axis(mol_idx=0,
                          atom_idx=demo.get_centroid_of_atoms(1, 2, mol_idx=0),
                          axis='x')
        for i in range(demo.mols[0].GetNumAtoms()):
            self.assertAlmostEqual(demo.get_point(i, 0).z, 0, -1)

    def pull_apart(self, mols: Chem.Mol, distance: float) -> Chem.Mol:
        walton = Walton(mols)
        walton.ring_on_plane(ring_idx=0, mol_idx=0)
        walton.align_by_mcs()
        walton.translate_parallel(mol_idx=1, distance=distance,
                                  base_atom_idx=0, pointer_atom_idx=2)
        walton(joining_cutoff=10)
        return walton.merged

    def test_pull_apart(self):
        benzene = Chem.MolFromSmiles('c1ccccc1')
        AllChem.EmbedMolecule(benzene)
        expectations = {0.0: 'c1ccccc1',  # benzene
                        2.5: 'c1ccc2ccccc2c1',  # naphthalene
                        3.5: 'C1CCC2(CC1)CCCCC2',  # spiro
                        4.0: 'c1ccc(-c2ccccc2)cc1', # bonded
                        4.5: 'c1ccc(Oc2ccccc2)cc1',  # oxy bridge
                        6.0: 'c1ccc(COc2ccccc2)cc1',  # methoxy bridge
                        7.0: 'c1ccc(CCOc2ccccc2)cc1',  # ethoxy bridge
                        8.0: 'c1ccc(OCCOc2ccccc2)cc1', # glycol bridge
                        }
        for distance, expected in expectations.items():
            merged: Chem.Mol = self.pull_apart([Chem.Mol(benzene), Chem.Mol(benzene)], distance)
            self.smiles_assertEqual(expected, merged)

    def test_tetrahedron(self):
        # build a tetrahedron!
        walton = Walton.from_smiles(ethanol='CCO', ethylamine='NCC')
        for i in range(2):
            walton.flatten_trio(mol_idx=i, atom_idcs=(0, 1, 2))
            walton.atom_to_origin(mol_idx=i, atom_idx=1)
        walton.rotate(mol_idx=1, theta=180, axis='z')
        walton.rotate(mol_idx=1, theta=90, axis='y')
        walton()
        self.smiles_assertEqual('CC(C)(N)O', walton.merged)

    def test_tetraherdron_with_atom_change(self):
        walton = Walton.from_smiles(methylether='COC', methylether2='COC')
        for i in range(2):
            walton.flatten_trio(mol_idx=i, atom_idcs=(0, 1, 2))
            walton.atom_to_origin(mol_idx=i, atom_idx=1)
        walton.rotate(mol_idx=i, theta=180, axis='z')
        walton.rotate(mol_idx=i, theta=90, axis='y')
        walton()
        self.smiles_assertEqual('CC(C)(C)C', walton.merged)
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

    def make_pair_by_split(self, conjoined: Chem.Mol, atom_idx: int) -> Tuple[Chem.Mol, Chem.Mol]:
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
        mol = Monster([benzene, moved]).combine().positioned_mol
        found = Chem.MolToSmiles(Chem.RemoveHs(mol))
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

from fragmenstein.monster.mcs_mapping import SpecialCompareAtoms

class Mappings(unittest.TestCase):
    """
    These test monster.mcs_mapping
    """
    toluene = Chem.MolFromSmiles('Cc1ccccc1')
    toluene.SetProp('_Name', 'toluene')

    benzyl = Chem.MolFromSmiles('*c1ccccc1')
    benzyl.SetProp('_Name', 'benzyl')

    methylpyrylium = Chem.MolFromSmiles('Cc1c[o+]ccc1')
    methylpyrylium.SetProp('_Name', 'methylpyrylium')
    # Monster.draw_nicely(None, methylpyrylium)

    methylpyridine = Chem.MolFromSmiles('Cc1ncccc1')
    methylpyridine.SetProp('_Name', 'methylpyridine')
    # Monster.draw_nicely(None, methylpyridine)


    def test_not_to_dummy(self):
        from rdkit.Chem import rdFMCS
        params = rdFMCS.MCSParameters()
        params.BondTyper = rdFMCS.BondCompare.CompareAny
        params.AtomTyper = SpecialCompareAtoms()
        compare = [self.benzyl, self.toluene]
        # hit -> followup
        res: rdFMCS.MCSResult = rdFMCS.FindMCS(compare, params)
        self.assertEqual(res.numAtoms, 6)  # there are 7 atoms, but only 6 are mapped as the dummy is excluded

    def test_user_map(self):
        """
        Map the azo with the oxo in pyridine and pyrylium.
        Using the methyl versions, the methyl does not map.
        """
        from rdkit.Chem import rdFMCS
        custom_map = {'methylpyridine': {2: 3}}
        params = rdFMCS.MCSParameters()
        params.BondTyper = rdFMCS.BondCompare.CompareAny
        params.AtomTyper = SpecialCompareAtoms(custom_map=custom_map)
        compare = [self.methylpyridine, self.methylpyrylium]
        # hit -> followup
        res: rdFMCS.MCSResult = rdFMCS.FindMCS(compare, params)
        self.assertEqual(res.numAtoms, 6, 'Mapping  the azo with the oxo means no methyl res:' +
                                         f' 6 matched, not {res.numAtoms}')
        full_maps = params.AtomTyper.get_valid_matches(params.AtomCompareParameters,
                                                      Chem.MolFromSmarts(res.smartsString),
                                                      self.methylpyridine, self.methylpyrylium
                                                      )
        # this is wrong: you can map the methyls to either side (ortho) of the heteroatom
        # self.assertEqual(len(full_maps), 1, f'There is only one way to map it not {len(full_maps)} ({full_maps})')
        for full_map in full_maps:
            for h, f in custom_map['methylpyridine'].items():
                self.assertEqual(dict(full_map)[h], f, f'hit atom {h} => {f} not {dict(full_map)[h]}')

    def test_user_negmap(self):
        """
        Ban index 1 on the azo=oxo mapping
        """
        from rdkit.Chem import rdFMCS
        custom_map = {'methylpyridine': {2: 3, 1: -2}}
        params = rdFMCS.MCSParameters()
        params.BondTyper = rdFMCS.BondCompare.CompareAny
        params.AtomTyper = SpecialCompareAtoms(custom_map=custom_map)
        compare = [self.methylpyridine, self.methylpyrylium]
        # hit -> followup
        res: rdFMCS.MCSResult = rdFMCS.FindMCS(compare, params)
        self.assertEqual(res.numAtoms, 5, 'ought to be azo=oxo and no index 1')
        full_maps = params.AtomTyper.get_valid_matches(params.AtomCompareParameters,
                                                       Chem.MolFromSmarts(res.smartsString),
                                                       self.methylpyridine, self.methylpyrylium
                                                       )
        # this is wrong: you can map the methyls to either side (ortho) of the heteroatom
        # self.assertEqual(len(full_maps), 1, f'There is only one way to map it not {len(full_maps)} ({full_maps})')
        for full_map in full_maps:
            self.assertNotIn(1, dict(full_map).keys())

    def test_flipper(self):
        from fragmenstein.monster.mcs_mapping import flip_mapping
        # Sequence[Tuple[int, int]]
        flipped = flip_mapping([(1, 2)])
        self.assertIsInstance(flipped, (list, tuple))


    def test_user_negmap2(self):
        """
        Assign followup index 1 to another molecule ('nullium' as its nothing) in the azo=oxo mapping
        thus preventing its mapping to methylpyridine
        """
        from rdkit.Chem import rdFMCS
        custom_map = {'methylpyridine': {2: 3}, 'nullium': {1: 1}}
        params = rdFMCS.MCSParameters()
        params.BondTyper = rdFMCS.BondCompare.CompareAny
        params.AtomTyper = SpecialCompareAtoms(custom_map=custom_map)
        compare = [self.methylpyridine, self.methylpyrylium]
        # hit -> followup
        res: rdFMCS.MCSResult = rdFMCS.FindMCS(compare, params)
        self.assertEqual(res.numAtoms, 5, 'ought to be azo=oxo and no index 1')
        full_maps = params.AtomTyper.get_valid_matches(params.AtomCompareParameters,
                                                       Chem.MolFromSmarts(res.smartsString),
                                                       self.methylpyridine, self.methylpyrylium
                                                       )
        # this is wrong: you can map the methyls to either side (ortho) of the heteroatom
        # self.assertEqual(len(full_maps), 1, f'There is only one way to map it not {len(full_maps)} ({full_maps})')
        for full_map in full_maps:
            self.assertNotIn(1, dict(full_map).values())


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
