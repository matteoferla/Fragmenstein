import logging
import tempfile
import unittest, os
# ======================================================================================================================
from multiprocessing import Process

import pyrosetta

pyrosetta.init(
    extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')
# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem

from fragmenstein import Monster, Victor, Igor, mpro_data, Walton
from fragmenstein.mpro import MProVictor
from typing import *
import numpy as np

# ======================================================================================================================


class MProPlaceTester(unittest.TestCase):

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
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneS = Chem.MolFromMolFile('test_mols/transtoluene2.mol')
        transtolueneS.SetProp('_Name', 'transtoluene-spiro')
        # cmd.rotate('z', -90, 'rototoluene', camera=0)  # this predates capt. Robert Walton
        mol = Monster(hits=[toluene, transtolueneS]).combine(keep_all=True).positioned_mol
        gotten = Chem.MolToSmiles(Chem.RemoveAllHs(mol))
        self.assertIn(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_real_merger(self):
        # Victor.enable_stdout(logging.DEBUG)
        x0138 = Chem.MolFromMolFile('test_mols/mac-x0138.mol')
        x0398 = Chem.MolFromMolFile('test_mols/mac-x0398.mol')
        monster = Monster([x0398, x0138])
        monster.combine()
        self.assertEqual('Nc1nc2c3c(c(O)cc(N)c3n1)C(O)=N2',
                         Chem.MolToSmiles(Chem.RemoveHs(monster.positioned_mol)),
                         )


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
        mol = Chem.MolFromMolFile('test_mols/F584.mol')
        flipped_F584 = 'COc1cccc2C(=O)NCCCc12'
        monster = Monster(hits=[hit_F584, ])
        flipped_F584 = 'COc1cccc2C(=O)NCCCc12'
        monster.place(Chem.MolFromSmiles(flipped_F584))
        # monster.show_comparison()
        # the problem is that a hydrogen gets mapped to a oxygen and this is unxpected (c.f. `_get_atom_maps`)
        self.assertEqual(len(monster.get_mcs_mappings(monster.initial_mol, monster.hits[0])[0][0]), 13)


class MultivictorPlaceTests(unittest.TestCase):
    def test_multivictor(self):
        from fragmenstein import MultiVictorPlacement
        to_place = Chem.MolFromMolFile('test_mols/placed_example1.mol')
        pdb_filename = 'test_mols/apo_example1.pdb'
        smiles = Chem.MolToSmiles(to_place)
        hits = [Chem.MolFromMolFile(os.path.join('test_mols', basename)) for basename in
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
            self.assertAlmostEqual(demo.get_point(i, 0).z, 0, 0)

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
