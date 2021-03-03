import unittest, os
from itertools import cycle, chain

# ======================================================================================================================
from typing import List

import pyrosetta
pyrosetta.init(
    extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')
# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem

from fragmenstein import Monster, Victor, Igor
from fragmenstein.mpro import MProVictor

# ======================================================================================================================
'''
To execute a single test:
python -m unittest tests.testDeLinker.DeLinkerTester.test_easy

To execute all test in a directory

python -m unittest discover tests/

'''

#TODO. Make this work 
class DeLinkerTester(unittest.TestCase):

    mPro_pdb = os.path.abspath( os.path.join(MProVictor.get_mpro_path(),  "template.pdb") )

    @classmethod
    def load_hits(cls, hit_codes):
      return  [ MProVictor.get_mol(xcode) for xcode in hit_codes ]

    def test_easy(self): #TODO. Breaks
        template = DeLinkerTester.mPro_pdb
        hit_codes= ['x0104', 'x1458']
        hits = self.load_hits(hit_codes= hit_codes )
        print( [ Chem.MolToSmiles(hit) for hit in hits])
        name = 'DeLinker_'+"-".join(hit_codes)
        Victor.monster_merging_mode="full"
        Victor.monster_joining_cutoff=5
        Victor.quick_renanimation = True
        victor = Victor(hits= hits, pdb_filename=template,
                           covalent_resi='145A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='CYS')
        victor.combine()
        gotten = Chem.MolToSmiles(Chem.RemoveHs(victor.minimised_mol))
        victor.make_pse()

    @classmethod
    def make_pse_generic(cls, mols: List[Chem.Mol], filename='test.pse', pdbFanme=None):
            """
            This is specifically for debugging the full fragment merging mode.
            For general use. Please use the Victor method ``make_pse``.
            :param mols: A list of molecules to display in pymol
            :param filename: pse file name
            :param tints: The colors for the molecules
            :return:
            """
            import pymol2
            tints = cycle(['wheat', 'palegreen', 'lightblue', 'paleyellow', 'lightpink', 'palecyan', 'lightorange',
                          'bluewhite'])
            assert '.pse' in filename, 'Must be a pymol pse extension!'
            with pymol2.PyMOL() as pymol:
                if pdbFanme:
                    pymol.cmd.load(pdbFanme)
                for i, mol in enumerate(mols):
                    try:
                        name = mol.GetProp('_Name')
                    except KeyError:
                        name= "mol_%d"%i
                    pymol.cmd.read_molstr(Chem.MolToMolBlock(mol, kekulize=False), name)
                    pymol.cmd.color(next(tints), f'{name} and name C*')
                pymol.cmd.save(filename)

    def test_no_linker_required(self):
        name = 'peridimethylnaphthalene'
        after = 'Cc1cccc2cccc(C)c12'
        template = os.path.join(MProVictor.get_mpro_path(), 'template.pdb')
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneF = Chem.MolFromMolFile('test_mols/transtoluene.mol')
        transtolueneF.SetProp('_Name', 'transtoluene-fuse')
        victor = Victor(hits=[toluene, transtolueneF],
                           pdb_filename=template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')
        victor.combine()
        gotten = Chem.MolToSmiles(Chem.RemoveHs(victor.minimised_mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')
        victor.make_pse()

    def test_pse(self):
      from joblib import load
      from random import sample
      generated_mols = load("/home/ruben/oxford/tools/Fragmenstein/fragmenstein/external/DeLinker_to_remove/molecules_generated.pkl")

      type(self).make_pse_generic(sample(list( chain.from_iterable(generated_mols[0].values())), 25),
                                  "/home/ruben/oxford/tools/Fragmenstein/fragmenstein/external/DeLinker_to_remove/molecules_generated.pse", pdbFanme=type(self).mPro_pdb)

    # def test_easy(self):
    #
    #     hits= self.load_hits(hit_codes=['x0107' ,'x0434', 'x1382'])
    #     Victor.quick_renanimation = True
    #     Victor.monster_merging_mode = "full"
    #     victor = Victor(  smiles = 'Cc1ccncc1NC(=O)Cc1cccc(Cl)c1',
    #                       hits = hits,
    #                       pdb_filename = type(self).mPro_pdb,
    #                       long_name = 'DeLinker_1',
    #                       ligand_resn = 'LIG',
    #                       covalent_resi='81A',
    #                       covalent_resn='CYS'
    #                       )
    #
    #     self.assertEqual(victor.error, '',victor.error)
    #     self.assertIsNotNone(victor.minimised_mol, 'Failed minimisation')
    #     victor.make_pse()

