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

# ======================================================================================================================
'''
To execute a single test:
python -m unittest test.MProTargetTester.test_easy

To execute all test in a directory

python -m unittest discover path/to/tests/dir

'''

class MProTargetTester(unittest.TestCase):

    def test_easy(self):
        MProVictor.quick_renanimation = True
        victor = MProVictor.from_hit_codes(smiles='CCNc1ncc(C#N)cc1CN1CCN(C(=O)C*)CC1',
                                  hit_codes=['x0692', 'x0305', 'x1249'],
                                  long_name='2_ACL')
        self.assertEqual(victor.error, '',victor.error)
        self.assertIsNotNone(victor.minimised_mol, 'Failed minimisation')
        victor.make_pse()

    def test_nasty(self):
        MProVictor.quick_renanimation = True
        victor = MProVictor.from_hit_codes(smiles='*CCC(=O)N1CC(CCN(C(=O)Nc2c(C)ncc(C)c2CCN2CCOCC2)c2cc(C)ccn2)C1',
                                           hit_codes=['x0434', 'x0540'],
                                           long_name='AGN-NEW-5f0-1_ACR1')
        self.assertEqual(victor.error, '', victor.error)
        self.assertIsNotNone(victor.minimised_mol, 'Failed minimisation')
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
                                           #hit_codes=['x0107','x0434','x0678','x0995','x1382'],
                                           hit_codes=['x0107' ,'x0434', 'x1382'],
                                           long_name='TRY-UNI-714a760b-6')
        self.assertEqual(victor.error, '', victor.error)
        self.assertIsNotNone(victor.minimised_mol, 'Failed minimisation')
        actual = MProVictor.get_mol('x2646')
        victor.make_pse(extra_mols=[actual])
        rmsd = victor.validate(reference_mol=actual)
        self.assertLess(rmsd, 1, f'The RMSD is large...')
        #self.assertIn('x1382', victor.monster.matched)
        #self.assertIn('x0995', victor.monster.unmatched) # red herring

# ======================================================================================================================


class RectifierTester(unittest.TestCase):
    def test_rectifier(self):
        # name: [before, after]
        chemdex = {'phenylnaphthalene': ('c1ccc2ccccc2c1(c3ccccc3)', 'c1ccc(-c2cccc3ccccc23)cc1'),
                   'benzo-azetine': ('C12CCCCC1CC2', 'C1CCC2CCCC2C1'),
                   #'conjoined': ('C1C2CCC2C1', 'C1CCCCC1'), # bridged hexane
                   'allene': ('C=C=C', 'C=CC'),
                   'benzo-cyclopronane': ('C12CCCCC1C2', 'C1CCC2CCCC2C1'),
                   #'norbornane': ('C1CC2CCC1C2', 'C1CC2CCC1C2'),
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
        mol.GetBondBetweenAtoms(0,1).SetBondType(Chem.BondType.SINGLE)
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
        mol.GetBondBetweenAtoms(0,1).SetBondType(Chem.BondType.SINGLE)
        mol.GetBondBetweenAtoms(6,7).SetBondType(Chem.BondType.AROMATIC)
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
        v = Victor.combine(hits=[toluene, rototoluene],
                           pdb_filename=template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')
        gotten = Chem.MolToSmiles(Chem.RemoveHs(v.minimised_mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_peridimethylnaphthalene(self):
        name = 'peridimethylnaphthalene'
        after = 'Cc1cccc2cccc(C)c12'
        template = os.path.join(MProVictor.get_mpro_path(), 'template.pdb')
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneF = Chem.MolFromMolFile('test_mols/transtoluene.mol')
        transtolueneF.SetProp('_Name', 'transtoluene-fuse')
        v = Victor.combine(hits=[toluene, transtolueneF],
                           pdb_filename=template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')
        gotten = Chem.MolToSmiles(Chem.RemoveHs(v.minimised_mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_spirodituluene(self):
        name = 'spirodituluene'
        after = 'C[C@@H]1C=CC[C@]2(C=C[C@H](C)C=C2)C1'
        template = os.path.join(MProVictor.get_mpro_path(), 'template.pdb')
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneS = Chem.MolFromMolFile('test_mols/transtoluene2.mol')
        transtolueneS.SetProp('_Name', 'transtoluene-spiro')
        # cmd.rotate('z', -90, 'rototoluene', camera=0)
        v = Victor.combine(hits=[toluene, transtolueneS],
                           pdb_filename=template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')
        gotten = Chem.MolToSmiles(Chem.RemoveHs(v.minimised_mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')
