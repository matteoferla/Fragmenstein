import logging
import unittest
import tempfile
from fragmenstein.demo import Mac1
from fragmenstein import Laboratory
import pandas as pd
import os

Laboratory.Victor.enable_stdout(logging.CRITICAL)

class LabTests(unittest.TestCase):
    def test_lab_combine(self):
        pdb_block = Mac1.get_template()
        hits = [Mac1.get_mol(f'diamond-{name}') for name in ['x0282_A','x0104_A','x0722_A','x0591_A','x0091_B']]
        lab = Laboratory(pdbblock=pdb_block, covalent_resi=None, run_plip=True)
        combinations: pd.DataFrame = lab.combine(hits, n_cores=8)
        self.assertGreater(len(combinations.outcome == 'acceptable'), 0)



