import unittest
import subprocess
import os

import pandas as pd

from fragmenstein import Victor
from fragmenstein.cli import FragmensteinParser
from fragmenstein import Laboratory
import fragmenstein.demo.test_mols as mol_folder
from fragmenstein.demo import Mac1
from pathlib import Path
from rdkit import Chem


# from tempfile import NamedTemporaryFile

class CliTests(unittest.TestCase):

    def command(self, command):
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        return out, err, process.returncode

    def test_monster_combine(self):
        parser = FragmensteinParser()
        return parser(['monster', 'combine',
                       '-i', str(Path(mol_folder.__file__).parent / 'mac-x0138.mol'),
                       str(Path(mol_folder.__file__).parent / 'mac-x0398.mol')
                       ])

    def test_monster_place(self):
        parser = FragmensteinParser()
        return parser(['monster', 'place',
                       '-i', str(Path(mol_folder.__file__).parent / 'mac-x0138.mol'),
                       str(Path(mol_folder.__file__).parent / 'mac-x0398.mol'),
                       '-s', 'COCCON'
                       ])

    def test_pipeline(self):
        # by default the results are written to the current working directory, yet we dont care here
        pdb_filename = (Path(Victor.work_path) / 'test.pdb').absolute().as_posix()
        sd_filename = (Path(Victor.work_path) / 'test.sdf').absolute().as_posix()
        os.chdir(Victor.work_path)
        print(  os.getcwd()  )
        # create the input files
        pdb_block = Mac1.get_template()
        hits = [Mac1.get_mol(f'diamond-{name}') for name in ['x0282_A', 'x0104_A', 'x0722_A']]  # , 'x0591_A', 'x0091_B'
        with Chem.SDWriter(sd_filename) as sdfh:
            for hit in hits:
                sdfh.write(hit)
        with open(pdb_filename, 'w') as pdbfh:
            pdbfh.write(pdb_block)
        # run the pipeline
        parser = FragmensteinParser()
        parser(['pipeline',
                '--input', sd_filename,
                '--template', pdb_filename,
                '--suffix', 'test',
                '--sw_length', '5',
                '--sw_databases', 'REAL-Database-22Q1.smi.anon'
                ])

    def zest_temp(self):
        """
        This was broken...
        """
        df_filename = (Path(Victor.work_path) / 'fragmenstein_placedtest.pkl.gz').absolute().as_posix()
        all_placements = pd.read_pickle(df_filename)
        Laboratory.export_sdf(df=all_placements)