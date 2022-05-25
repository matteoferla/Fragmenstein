import unittest

import numpy as np
# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem

from fragmenstein import Monster, Victor, Igor
from fragmenstein.demo import TestSet, MPro


# ======================================================================================================================

# ======================================================================================================================

class MonsterPlaceTests(unittest.TestCase):

    def test_sample_new_conformation(self):
        smiles = "C1C2C(C=C(C=2)C(C2C=CC=C2)CNOC)C=CC=1"
        mol = Chem.MolFromSmiles(smiles)
        AllChem.EmbedMolecule(mol, randomSeed=131)
        frag_mol = Chem.FragmentOnBonds(mol, [5, 11], addDummies=False)
        hits = Chem.GetMolFrags(frag_mol, asMols=True)[:2]
        # from matplotlib import pyplot as plt
        # from rdkit.Chem import Draw
        # plt.imshow(Draw.MolsToGridImage(hits)); plt.show()
        ori_monster = Monster(hits=hits, random_seed=131)
        ori_monster.place_smiles(smiles)
        # sample_new_conformation
        seeds = [121, 23421, 1]
        for se1 in seeds:
            mol = ori_monster.sample_new_conformation(random_seed=se1)
            coords1 = mol.GetConformer().GetPositions()
            for se2 in seeds:
                mol = ori_monster.sample_new_conformation(random_seed=se2)
                coords2 = mol.GetConformer().GetPositions()
                if se1 == se2:
                    self.assertAlmostEqual(np.sum(np.abs(coords1 - coords2)), 0)
                else:
                    self.assertTrue(np.sum(np.abs(coords1 - coords2)) > 3)

    def test_random_seed(self):
        smiles = "C1C2C(C=C(C=2)C(C2C=CC=C2)CNOC)C=CC=1"
        mol = Chem.MolFromSmiles(smiles)
        AllChem.EmbedMolecule(mol, randomSeed=131)
        frag_mol = Chem.FragmentOnBonds(mol, [5, 11], addDummies=False)
        hits = Chem.GetMolFrags(frag_mol, asMols=True)[:2]
        # from matplotlib import pyplot as plt
        # from rdkit.Chem import Draw
        # plt.imshow(Draw.MolsToGridImage(hits)); plt.show()
        seeds = [121, 23421, 1]
        for se1 in seeds:
            mol = Monster(hits=hits, random_seed=se1).place_smiles(smiles).positioned_mol
            coords1 = mol.GetConformer().GetPositions()
            for se2 in seeds:
                mol = Monster(hits=hits, random_seed=se2).place_smiles(smiles).positioned_mol
                coords2 = mol.GetConformer().GetPositions()
                if se1 == se2:
                    self.assertAlmostEqual(np.sum(np.abs(coords1 - coords2)), 0)
                else:
                    self.assertTrue(np.sum(np.abs(coords1 - coords2)) > 3)

    def test_flipped_lactam(self):
        """
        Given a benzo + 7-membered lactam map a mol with the amide flipped
        """
        hit_F584 = TestSet.get_mol('F584')
        monster = Monster(hits=[hit_F584, ])
        monster.place_smiles(smiles='COc1cccc2C(=O)NCCCc12', long_name='flipped_F584')
        # monster.show_comparison()
        self.assertEqual(len(monster.get_mcs_mappings(followup=monster.initial_mol, hit=monster.hits[0])[0][0]), 13)

    def test_forced_flipped_lactam(self):

        """
        Given a benzo + 7-membered lactam map a mol with the amide flipped as the previous test
        force the ketones to match
        """
        hit_F584 = TestSet.get_mol('F584')
        hit_F584.SetProp('_Name', 'hit_F584')
        monster = Monster(hits=[hit_F584, ])
        monster.place_smiles(smiles='COc1cccc2C(=O)NCCCc12',
                             long_name='flipped_F584',
                             custom_map={'hit_F584': {13: 8}},
                             # merging_mode='expansion', # works on both expansion and no blend modes
                             )
        # monster.draw_nicely(monster.hits[0])
        # monster.draw_nicely(monster.initial_mol)
        # monster.show_comparison()
        # the [0][0][13] is because of:
        # Tuple[List[Dict[int, int]], ExtendedFMCSMode]
        self.assertIn(13, monster.get_mcs_mappings(followup=monster.initial_mol, hit=monster.hits[0])[0][0])
        self.assertEqual(monster.get_mcs_mappings(followup=monster.initial_mol, hit=monster.hits[0])[0][0][13], 8)
        # self.assertEqual(len(monster.get_mcs_mappings(followup=monster.initial_mol, hit=monster.hits[0])[0][0]), 13)

    def test_thiazole_flip(self):
        """
        See ``.get_5SB7_mols`` for details.

        Potentially the number of specified atoms is excessive,
        but as a human I just wrote down the atoms in the ring and was done with it.
        """
        mols = TestSet.get_5SB7_mols()
        monster = Monster([mols[0]])
        monster.place(Chem.Mol(mols[0]),
                      merging_mode='expansion',
                      # custom_map={'F36': {1:7}, 'F04': {4:7}}
                      custom_map={'F04': {-1: 4,  # no amine
                                          4: -2,  # no amine
                                          12: 12,
                                          6: 13,
                                          13: 6,
                                          15: 14,
                                          14: 15}}
                      )
        self.assertEqual(len(list(filter(len, monster.origin_from_mol(monster.positioned_mol)))), 15)
        # the amine is banned in the map and the ring is flipped.

    def test_thiazole_followup(self):
        """
        The followup compound is from 5SB3
        and is 'bold4' in mols.
        """
        mols = TestSet.get_5SB7_mols()
        monster = Monster(mols[:2])
        monster.place(Chem.Mol(mols[3]),
                      merging_mode='expansion',
                      custom_map={'F36': {1: 7},
                                  'F04': {4: -1,  # no amine
                                          12: 13,  # root to Ph
                                          13: 5,
                                          6: 14,
                                          }
                                  }
                      )
        # monster.show_comparison()
        # monster.to_nglview(True)
        # the amine at 7 is banned in F04.
        self.assertEqual(monster.origin_from_mol()[7], ['F36.1'])
        self.assertEqual(len(list(filter(len, monster.origin_from_mol(monster.positioned_mol)))), 22)

    def test_renumber(self):
        mols = TestSet.get_5SB7_mols()
        # F04 {root: 12, N: 13, amine: 0, C: 6},
        # original {root: 13, N: 5, amine: 7, C: 14},
        # renumbered {root: 8, N:19, amine: 12, C: 9}
        renumbered = Monster.renumber_followup_custom_map(mols[3],
                                                          Chem.MolFromSmiles(Chem.MolToSmiles(mols[3])),
                                                          custom_map={'F36': {1: 7},
                                                                      'F04': {0: -1,  # no amine
                                                                              -1: 7,
                                                                              13: 5,  # to N
                                                                              12: 13,  # root to Ph
                                                                              6: 14,  # to C
                                                                              }
                                                                      }

                                                          )
        expected = {'F36': {1: 12}, 'F04': {0: -1, -1: 12, 13: 19, 12: 8, 6: 9}}
        for name in renumbered:
            self.assertEqual(tuple(renumbered[name]), tuple(expected[name]),
                             f'{name} {renumbered[name]} != {expected[name]}')

    def test_victor_mol(self):
        """
        Victor from mol has different indices than from a smiles
        """
        mols = TestSet.get_5SB7_mols()
        Igor.init_pyrosetta()

        victor = Victor(hits=mols[:2], pdb_block=MPro.get_template(), ligand_resi='1X')

        victor.place(mols[3],
                     long_name=mols[3].GetProp('_Name'),
                     merging_mode='expansion',
                     custom_map={'F36': {1: 7},
                                 'F04': {0: -1,  # no amine
                                         -1: 7,  # no amine
                                         13: 5,  # to N
                                         12: 13,  # root to Ph
                                         6: 14,  # to C
                                         }
                                 }
                     )
        # victor.show_comparison()
        # victor.to_nglview()
        rsmd = victor.validate(mols[3])['reference2minimized_rmsd']
        self.assertLess(rsmd, 1.1, f"The resulting RMSD from the crystal is {rsmd}, which is greater than 1.")


if __name__ == '__main__':
    unittest.main()
