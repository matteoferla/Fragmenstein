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

from fragmenstein.demo import TestSet

# import logging
# Victor.enable_stdout(logging.DEBUG)
# ======================================================================================================================

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
        toluene = TestSet.get_mol('toluene')
        rototoluene = TestSet.get_mol('rototoluene')
        mol = Monster(hits=[toluene, rototoluene]).combine(keep_all=True).positioned_mol
        gotten = Chem.MolToSmiles(Chem.RemoveHs(mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_peridimethylnaphthalene(self):
        name = 'peridimethylnaphthalene'
        after = 'Cc1cccc2cccc(C)c12'
        toluene = TestSet.get_mol('toluene')
        transtolueneF = TestSet.get_mol('transtoluene')
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
        toluene = TestSet.get_mol('toluene')
        transtolueneS = TestSet.get_mol('transtoluene2')
        transtolueneS.SetProp('_Name', 'transtoluene-spiro')
        # cmd.rotate('z', -90, 'rototoluene', camera=0)  # this predates capt. Robert Walton
        mol = Monster(hits=[toluene, transtolueneS]).combine(keep_all=True).positioned_mol
        gotten = Chem.MolToSmiles(Chem.RemoveAllHs(mol))
        self.assertIn(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_real_merger(self):
        # Victor.enable_stdout(logging.DEBUG)
        x0138 = TestSet.get_mol('mac-x0138')
        x0398 = TestSet.get_mol('mac-x0398')
        monster = Monster([x0398, x0138])
        monster.combine()
        self.assertEqual('Nc1nc2c3c(c(O)cc(N)c3n1)C(O)=N2',
                         Chem.MolToSmiles(Chem.RemoveHs(monster.positioned_mol)),
                         )
