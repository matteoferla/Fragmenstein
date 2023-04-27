import unittest
from typing import *

from fragmenstein import Victor, Igor, mpro_data
from fragmenstein.mpro import MProVictor


# ======================================================================================================================

class MProPlaceTester(unittest.TestCase):

    def setUp(self):
        Igor.init_pyrosetta()

    def untest_red_herring(self):  # without the test_ word this will not run.
        """
        To a **human** this looks easy. x0692 is a red herring and the other two make the molecule.
        As currently written this will fail.

        See red herring test notes in documentation.

        :return:
        """
        # PAU-UNI-52c0427f-1
        MProVictor.quick_reanimation = True
        victor = MProVictor.from_hit_codes(hit_codes=['x0692', 'x0305', 'x1249'])
        victor.place('CCNc1ncc(C#N)cc1CN1CCN(C(=O)C*)CC1', long_name='2_ACL')
        self.assertEqual(victor.error_msg, '', victor.error_msg)
        self.assertIsNotNone(victor.minimized_mol, 'Failed minimisation')
        msg = f'x1249 is the red herring, prediction: {victor.monster.unmatched} discarded, ' + \
              f'while x0305 and x0692 are the true inspirations. kept: {victor.monster.matched}'
        self.assertIn('x1249', victor.monster.unmatched, msg)
        self.assertIn('x0305', victor.monster.matched, msg)

        victor.make_pse()

    def untest_nasty(self):  # without the test_ word this will not run.
        """
        The human suggested a lot of novel groups.
        'x0540' is a really odd inspiration. Three atoms are conserved. the rest aren't.

        Like the red herring this is impossible.
        :return:
        """
        MProVictor.quick_reanimation = True
        victor = MProVictor.from_hit_codes(hit_codes=['x0434', 'x0540'])
        victor.place('*CCC(=O)N1CC(CCN(C(=O)Nc2c(C)ncc(C)c2CCN2CCOCC2)c2cc(C)ccn2)C1',
                     long_name='AGN-NEW-5f0-1_ACR1')
        self.assertEqual(str(victor.error_msg), '', str(victor.error_msg))
        self.assertIsNotNone(victor.minimized_mol, 'Failed minimisation')
        self.assertEqual(len(victor.monster.unmatched), 0,
                         f'Both were correct but {victor.monster.unmatched} was discarded')
        victor.make_pse()

    def test_incorrect(self):
        """
        This case has two hits that are not ideal. One, x0995, totally barney.
        These will be rejected.

        :return:
        """
        MProVictor.quick_reanimation = True
        victor = MProVictor.from_hit_codes(hit_codes='x0692,x0770,x0995'.split(','))
        victor.place('*C(=N)CN1CCN(Cc2ccc(-c3cc(CC)ncn3)c(F)c2)CC1',
                     long_name='BEN-VAN-c98-4')
        self.assertEqual(victor.error_msg, '', victor.error_msg)
        self.assertIsNotNone(victor.minimized_mol, 'Failed minimisation')
        # in unmatched, x0995 is the red herring.
        self.assertIn('x0995', victor.monster.unmatched, f'x0995 is the red herring')
        #victor.make_pse()

    def test_pentachromatic(self):
        """
        This hit fails to identify that the extra chloride comes from x1382.
        The human chose these.

        * the methylpyridine off x0107
        * the benzene-pyridine off x0434, but wanted an amide not a ureido
        * the link between the rings as amide x0678
        * x0995 is a red herring
        * benzene with a chroride from x1382
        """
        MProVictor.quick_reanimation = True
        # ,'x2646'
        Victor.monster_throw_on_discard = True
        victor = MProVictor.from_hit_codes(  # hit_codes=['x0107','x0434','x0678','x0995','x1382'],
            hit_codes=['x0107', 'x0434', 'x1382'])
        victor.place('Cc1ccncc1NC(=O)Cc1cccc(Cl)c1',
                     long_name='TRY-UNI-714a760b-6')
        self.assertEqual(victor.error_msg, '', victor.error_msg)
        self.assertIsNotNone(victor.minimized_mol, 'Failed minimisation')
        actual = mpro_data.get_mol('x2646')
        try:
            victor.make_pse(extra_mols=[actual])
        except ModuleNotFoundError:
            print('PyMOL is not installed. Cannot export result.')
        validation: Dict[str, float] = victor.validate(reference_mol=actual)
        rmsd = validation['reference2minimized_rmsd']
        self.assertLess(rmsd, 2, f'The RMSD is large...')
        self.assertIn('x1382', victor.monster.matched)
        self.assertNotIn('x0995', victor.monster.unmatched)  # red herring

if __name__ == '__main__':
    unittest.main()
