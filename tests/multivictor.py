import os
import tempfile
import unittest

# ======================================================================================================================
from rdkit import Chem

from fragmenstein import Victor, Igor
from fragmenstein.demo import TestSet

# ======================================================================================================================

class MultivictorPlaceTests(unittest.TestCase):

    def setUp(self):
        Igor.init_pyrosetta()

    def test_multivictor(self):
        from fragmenstein import MultiVictorPlacement
        to_place = TestSet.get_mol('placed_example1')
        pdb_block = TestSet.get_text('apo_example1.pdb')
        smiles = Chem.MolToSmiles(to_place)
        hits = [TestSet.get_mol('x0032_0A'),
                #TestSet.get_mol('x0103_0A')
                ]
        # Victor.enable_stdout(level=logging.ERROR)
        with tempfile.TemporaryDirectory() as tmpdir:
            Victor.work_path = os.path.join(tmpdir, "multivictor_out")
            mv = MultiVictorPlacement(hits=hits, pdb_block=pdb_block)
            mv.place(smiles, number_runs=4)
            # print(mv.retrieve_best_victor())
            # print(mv.retrieve_scores())
            self.assertLess(mv.retrieve_scores()[0], -7)


if __name__ == '__main__':
    unittest.main()
