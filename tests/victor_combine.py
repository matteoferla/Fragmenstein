import unittest
from fragmenstein import Victor, Igor
from fragmenstein.demo import MPro, Mac1
from rdkit import Chem

class VictorCombineTests(unittest.TestCase):

    def setUp(self):
        Igor.init_pyrosetta()

    def test_merger(self):
        """
        This is a case wherein there is a 1B residue already.
        """
        hits = [MPro.get_mol('x0395'), MPro.get_mol('x0434')]
        vicky = Victor(hits, pdb_block=MPro.get_template())
        vicky.combine()
        self.assertEquals(len(vicky.summarize()['disregarded']), 0)

    def test_glitch_from_tim(self):
        hit_block1 = 'DLS-EU0046_0A\n     RDKit          3D\n\n  9  9  0  0  0  0  0  0  0  0999 V2000\n   30.3230   48.4150   -5.7320 N   0  0  0  0  0  0  0  0  0  0  0  0\n   29.7040   48.9230   -6.5500 C   0  0  0  0  0  0  0  0  0  0  0  0\n   28.9400   49.5100   -7.6020 C   0  0  0  0  0  0  0  0  0  0  0  0\n   28.2880   48.9010   -8.7080 C   0  0  0  0  0  0  0  0  0  0  0  0\n   28.1240   47.4460   -9.0030 C   0  0  0  0  0  0  0  0  0  0  0  0\n   28.7750   50.9000   -7.8120 C   0  0  0  0  0  0  0  0  0  0  0  0\n   29.2110   51.8940   -7.0630 N   0  0  0  0  0  0  0  0  0  0  0  0\n   28.0730   51.0510   -8.9820 N   0  0  0  0  0  0  0  0  0  0  0  0\n   27.7850   49.8260   -9.5300 N   0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  3  0\n  3  2  1  0\n  4  3  1  0\n  5  4  1  0\n  6  3  2  0\n  7  6  1  0\n  8  6  1  0\n  9  4  2  0\n  9  8  1  0\nM  END\n'
        hit_block2 = 'DLS-EU0056_0A\n     RDKit          3D\n\n  9  9  0  0  0  0  0  0  0  0999 V2000\n   28.3400   50.9740   -8.0450 Br  0  0  0  0  0  0  0  0  0  0  0  0\n   28.5950   52.7780   -7.5400 C   0  0  0  0  0  0  0  0  0  0  0  0\n   28.1520   53.7870   -8.3680 C   0  0  0  0  0  0  0  0  0  0  0  0\n   29.1390   53.0770   -6.3090 C   0  0  0  0  0  0  0  0  0  0  0  0\n   29.2440   54.3360   -5.8740 N   0  0  0  0  0  0  0  0  0  0  0  0\n   28.7990   55.3080   -6.6760 C   0  0  0  0  0  0  0  0  0  0  0  0\n   28.2350   55.1020   -7.9250 C   0  0  0  0  0  0  0  0  0  0  0  0\n   27.6630   56.2440   -8.7350 C   0  0  0  0  0  0  0  0  0  0  0  0\n   28.7420   56.9200   -9.4540 N   0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  0\n  3  2  2  0\n  4  2  1  0\n  5  4  2  0\n  6  5  1  0\n  7  3  1  0\n  7  6  2  0\n  8  7  1  0\n  9  8  1  0\nM  END\n'
        hits = [Chem.MolFromMolBlock(hit_block1), Chem.MolFromMolBlock(hit_block2)]
        vicky = Victor(hits, pdb_block=Mac1.get_template())
        vicky.combine()
        mol = vicky.minimized_mol
        self.assertFalse(bool(Chem.DetectChemistryProblems(mol)))
        print(Chem.MolToSmiles(mol))
        # bad: [H]c1nc([H])c(C([H])([H])N([H])[H])c([H])c1N([H])C1:N:N:C(C([H])([H])[H]):C:1C#N
        self.assertIsNotNone(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))


