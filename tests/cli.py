import unittest
import subprocess

from fragmenstein.cli import FragmensteinParser
import pkg_resources
from fragmenstein.demo import TestSet
import fragmenstein.demo.test_mols as mol_folder
from pathlib import Path

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





