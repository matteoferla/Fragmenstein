import logging
import tempfile
import unittest, os
# ======================================================================================================================
from multiprocessing import Process

from fragmenstein.igor.pyrosetta_import import pyrosetta

# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem, rdqueries
from rdkit.Geometry import Point3D

# TESTS IS EXTERNAL TO FRAGMENSTEIN DO NOT CHANGE TO RELATIVE!
from fragmenstein import Monster, Igor, Victor
from fragmenstein.demo import Mac1
from fragmenstein.error import FragmensteinError
from fragmenstein.mpro import MProVictor
from typing import *
import numpy as np
import json
from rdkit import Chem, rdBase
from fragmenstein import Laboratory, Victor

# ======================================================================================================================

class Internals(unittest.TestCase):

    def setUp(self):
        Victor.capture_rdkit_log()
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
        dummies = mol.GetAtomsMatchingQuery(rdqueries.AtomNumEqualsQueryAtom(0))
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
        except FragmensteinError as error:
            pass

    def test_neigh_bonding(self):
        """
        This was to address the bug from 5/1/23.

        Issue was in ``_nan_fill_submatrix``, which was not blanking due to
        ``isinstance(np.int64(42), int) == False``
        """
        # make mols
        methane: Chem.Mol = Chem.MolFromSmiles('C')
        methane.SetProp('_Name', 'methane')
        AllChem.EmbedMolecule(methane, coordMap={0: Point3D(0, 0, 0)})
        ammonia = Chem.MolFromSmiles('N')
        ammonia.SetProp('_Name', 'ammonia')
        # merge
        monster = Monster([])
        AllChem.EmbedMolecule(ammonia, coordMap={0: Point3D(3, 0, 0)})
        self.translate(ammonia, x=3)
        # monster.join_neighboring_mols(methane, ammonia)
        combo, candidates = monster._find_all_closest(methane, ammonia)  # _find_all_closest is in communal
        self.assertEqual(combo.GetNumAtoms(), 2)
        self.assertNotEqual(candidates[0][0], candidates[0][1])

    def test_percent_hybrid(self):
        chlorobutane = Chem.MolFromSmiles('[Cl]CCCC')
        origins = [[], ['ethanol.1'], ['ethanol.2', 'isopronanol.3'], ['isopronanol.2'], ['isopronanol.4']]
        chlorobutane.SetProp('_Origins', json.dumps(origins))
        hybrid = Laboratory.percent_hybrid(None, chlorobutane)
        self.assertEqual(hybrid, 34)  # ethanol has 1 single origin atom, isopronanol has 2

    def test_assign_empty(self):
        vicky = Victor.__new__(Victor)
        vicky.apo_pdbblock = Mac1.get_template()
        print(vicky._get_empty_resi())



    # def test_doubleconstraint(self):
    #     diaminopentane = Chem.MolFromSmiles('NCCCCCN')
    #     AllChem.EmbedMolecule(diaminopentane)
    #     AllChem
    #
    #     #Monster([])

if __name__ == '__main__':
    unittest.main()
