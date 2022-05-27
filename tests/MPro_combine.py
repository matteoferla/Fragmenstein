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



if __name__ == '__main__':
    unittest.main()
