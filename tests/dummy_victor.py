import logging
import tempfile
import unittest, os
# ======================================================================================================================
from multiprocessing import Process

import pyrosetta

# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem
from fragmenstein.demo import TestSet
from fragmenstein import Monster, Victor, Igor, mpro_data, Walton, demo
from typing import *
import numpy as np

# ======================================================================================================================

class DummyCombineTests(unittest.TestCase):
    def setUp(self):
        Igor.init_pyrosetta()

    def making_of_dummy(self):
        # This is the code I used to make the dummy.pdb file
        import numpy as np
        import pyrosetta
        from types import ModuleType
        prp: ModuleType = pyrosetta.rosetta.protocols
        prc: ModuleType = pyrosetta.rosetta.core
        pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
        import pyrosetta_help as ph

        Igor.init_pyrosetta()
        pose = pyrosetta.pose_from_sequence('ELVIS')
        scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
        prp.relax.FastRelax(scorefxn, 5).apply(pose)
        pose.translate(np.array([-100, 0, 0]))
        block = ph.get_pdbstr(pose)
        return block
    def test_ring_on_line(self):
        """
        This is a case to test what happens to a ring placed on a line:
        It does not open but it does not place with a low RMSD.
        """
        block: str = TestSet.get_text('dummy.pdb')
        hydroaminohexene = Chem.MolFromSmiles('OC=C/C=C/C=CN')
        hydroaminohexene.SetProp('_Name', 'hexane')
        AllChem.EmbedMolecule(hydroaminohexene)
        victor = Victor([hydroaminohexene], pdb_block=block)
        victor.place(Chem.MolFromSmiles('Oc1ccccc1N'))
        self.assertGreater(victor.summarize()['comRMSD'], 1, 'Ring on line did not shift')

if __name__ == '__main__':
    unittest.main()