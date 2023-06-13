import unittest

from fragmenstein import Igor
from fragmenstein.mpro import MProVictor


# ======================================================================================================================
# ======================================================================================================================


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

    def test_drift_prevention_glitch(self):
        """
        I am not sure why but I had to specify which atoms mapped to which in the drift prevention step
        else it gives ``RuntimeError: No sub-structure match found between the probe and query mol``
        """
        MProVictor.quick_reanimation = False
        victor = MProVictor.from_hit_codes(hit_codes=['x0426', 'x0540'])
        victor.combine()


if __name__ == '__main__':
    unittest.main()
