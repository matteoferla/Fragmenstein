import unittest, os
# ======================================================================================================================
import pyrosetta
pyrosetta.init(
    extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')
# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem

from fragmenstein import Monster, Victor, Igor, Rectifier
from fragmenstein.mpro import MProVictor

# ======================================================================================================================
'''
To execute a single test:
python -m unittest tests.testDeLinker.DeLinkerTester.test_easy

To execute all test in a directory

python -m unittest discover tests/

'''

class DeLinkerTester(unittest.TestCase):

    mPro_pdb = os.path.abspath( os.path.join(MProVictor.get_mpro_path(),  "template.pdb") )

    @classmethod
    def load_hits(cls, hit_codes):
      return  [ MProVictor.get_mol(xcode) for xcode in hit_codes ]

    def test_easy(self):
        template = os.path.join(MProVictor.get_mpro_path(), 'template.pdb')
        hit_codes= ['x0104', 'x1458']
        hits = self.load_hits(hit_codes= hit_codes )
        print( [ Chem.MolToSmiles(hit) for hit in hits])
        name = 'DeLinker-'+"-".join(hit_codes)
        Victor.monster_merging_mode="full"
        Victor.monster_joining_cutoff=10
        victor = Victor.combine(hits= hits,
                           pdb_filename=template,
                           covalent_resi='145A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='CYS')
        gotten = Chem.MolToSmiles(Chem.RemoveHs(victor.minimised_mol))
        victor.make_pse()


    def test_no_linker_required(self):
        name = 'peridimethylnaphthalene'
        after = 'Cc1cccc2cccc(C)c12'
        template = os.path.join(MProVictor.get_mpro_path(), 'template.pdb')
        toluene = Chem.MolFromMolFile('test_mols/toluene.mol')
        toluene.SetProp('_Name', 'toluene')
        transtolueneF = Chem.MolFromMolFile('test_mols/transtoluene.mol')
        transtolueneF.SetProp('_Name', 'transtoluene-fuse')
        victor = Victor.combine(hits=[toluene, transtolueneF],
                           pdb_filename=template,
                           covalent_resi='3A',  # a random residue is still required for the constaint ref atom.
                           covalent_resn='VAL')
        gotten = Chem.MolToSmiles(Chem.RemoveHs(victor.minimised_mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')
        victor.make_pse()

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

