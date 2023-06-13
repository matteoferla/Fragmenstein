# Testing discarding
from rdkit import Chem
from fragmenstein import Walton, Victor, Igor
import unittest
from fragmenstein import DistanceError

class TestDiscard(unittest.TestCase):
    pdb_block = \
        '''ATOM      1  N   ALA A   1      -0.677  -1.230  -0.491  1.00  0.00           N  
        ATOM      2  CA  ALA A   1      -0.001   0.064  -0.491  1.00  0.00           C  
        ATOM      3  C   ALA A   1       1.499  -0.110  -0.491  1.00  0.00           C  
        ATOM      4  O   ALA A   1       2.030  -1.227  -0.502  1.00  0.00           O  
        ATOM      5  CB  ALA A   1      -0.509   0.856   0.727  1.00  0.00           C  
        ATOM      6  H   ALA A   1      -0.131  -2.162  -0.491  1.00  0.00           H  
        ATOM      7  HA  ALA A   1      -0.269   0.603  -1.418  1.00  0.00           H  
        ATOM      8 1HB  ALA A   1      -1.605   1.006   0.691  1.00  0.00           H  
        ATOM      9 2HB  ALA A   1      -0.285   0.342   1.681  1.00  0.00           H  
        ATOM     10 3HB  ALA A   1      -0.053   1.861   0.784  1.00  0.00           H  
        TER   
        END'''

    def setUp(self) -> None:
        Igor.init_pyrosetta()

    def get_separated(self) -> Victor:
        wally = Walton.from_smiles(benzene='c1ccccc1', tetrazole='[nH]1nncc1')
        wally.ring_on_plane(ring_idx=0, mol_idx=0)
        wally.superpose_by_mcs()
        wally.translate_parallel(mol_idx=1,
                                 distance=10,
                                 base_atom_idx=0,
                                 pointer_atom_idx=2)
        Victor.monster_throw_on_discard = True
        vicky = Victor(wally.mols,
                       pdb_block=self.pdb_block,
                       ligand_resn='LIG')
        return vicky
    def test_discard(self):
        vicky: Victor = self.get_separated()
        try:
            vicky.combine(joining_cutoff=3)
            print(vicky.summarize())
            self.fail('Should have thrown')
        except DistanceError:
            print('error thrown as expected')

    def test_success(self):
        vicky: Victor = self.get_separated()
        vicky.monster.throw_on_discard = True
        vicky.combine(joining_cutoff=12)
        print(vicky.summarize())

    def test_success2(self):
        vicky: Victor = self.get_separated()
        vicky.place('c1nnn(CCCCCOc2ccccc2)c1')
        print(vicky.summarize())

# these test need a bit more meat:
# the second benzene is not a contributor as all overlap.
# wally = Walton.from_smiles(benzene='c1ccccc1', benzene2='c1ccccc1')
# wally.ring_on_plane(ring_idx=0, mol_idx=0)
# wally.superpose_by_mcs()
# vicky = Victor(wally.mols,
#                pdb_block=pdb_block,
#                ligand_resn='LIG')
# # vicky.monster_throw_on_discard = True <-- NO!
# vicky.monster.throw_on_discard = True
# vicky.place('c1cccc(CCCCCOc2ccccc2)c1', long_name='discard-error2')

