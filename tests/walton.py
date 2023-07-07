import unittest

# ======================================================================================================================
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List
from fragmenstein import Walton

# ======================================================================================================================

class WaltonTests(unittest.TestCase):

    def smiles_assertEqual(self, a, b):
        """
        helper method to test equality of smiles.
        Simply reshuffled the parts. Lazy.
        """
        if isinstance(a, str):
            expected_smiles = a
            mol = b
        else:
            expected_smiles = b
            mol = a
        obtained_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
        self.assertEqual(expected_smiles, obtained_smiles)

    def test_superpose_map(self):
        demo = Walton.from_smiles(resorcinol='c1ccc(O)cc1O', eugenol='Oc1ccc(cc1OC)CC=C')  # create instance
        demo.superpose_by_map({(0, 1): {4: 0, 3: 1, 2: 2}})  # superpose molecules by atom indices
        demo()  # merge (Fragmenstein's Monster)
        self.smiles_assertEqual('C=CCc1ccc(O)c(OC)c1O', demo.merged)

    def test_polygon(self):
        benzene = Walton.create_polygon(n=6, bond_type=Chem.BondType.AROMATIC)
        self.smiles_assertEqual('c1ccccc1', benzene)

    def test_ring_on_plane(self):
        demo = Walton.from_smiles(furan='c1ccco1')
        demo.ring_on_plane(mol_idx=0)  # flat on xy
        demo.atom_on_axis(mol_idx=0, atom_idx=4, axis='x')  # ox on axis
        # the atom_idx can be actually be passed a Geometry.Point3D:
        demo.atom_on_axis(mol_idx=0,
                          atom_idx=demo.get_centroid_of_atoms(1, 2, mol_idx=0),
                          axis='x')
        for i in range(demo.mols[0].GetNumAtoms()):
            self.assertAlmostEqual(demo.get_point(i, 0).z, 0, -1)

    def pull_apart(self, mols: List[Chem.Mol], distance: float) -> Chem.Mol:
        walton = Walton(mols)
        walton.ring_on_plane(ring_idx=0, mol_idx=0)
        walton.superpose_by_mcs()
        walton.translate_parallel(mol_idx=1, distance=distance,
                                  base_atom_idx=0, pointer_atom_idx=2)
        walton(joining_cutoff=10)
        return walton.merged

    def test_pull_apart(self):
        """
        Two benzene molecules are placed at different distances.
        """
        benzene = Chem.MolFromSmiles('c1ccccc1')
        benzene.SetProp('_Name', 'benzene')
        AllChem.EmbedMolecule(benzene)
        expectations = {0.0: 'c1ccccc1',  # benzene
                        2.5: 'c1ccc2ccccc2c1',  # naphthalene
                        #3.5: 'C1CCC2(CC1)CCCCC2',  # spiro... inconsistent!
                        4.0: 'c1ccc(-c2ccccc2)cc1', # bonded
                        4.5: 'c1ccc(Oc2ccccc2)cc1',  # oxy bridge
                        6.0: 'c1ccc(COc2ccccc2)cc1',  # methoxy bridge
                        7.0: 'c1ccc(CCOc2ccccc2)cc1',  # ethoxy bridge
                        8.0: 'c1ccc(OCCOc2ccccc2)cc1', # glycol bridge
                        }
        for distance, expected in expectations.items():
            print(distance)
            merged: Chem.Mol = self.pull_apart([Chem.Mol(benzene), Chem.Mol(benzene)], distance)
            self.smiles_assertEqual(expected, merged)

    def test_tetrahedron(self):
        # build a tetrahedron!
        walton = Walton.from_smiles(ethanol='CCO', ethylamine='NCC')
        for i in range(2):
            walton.flatten_trio(mol_idx=i, atom_idcs=(0, 1, 2))
            walton.atom_to_origin(mol_idx=i, atom_idx=1)
        walton.rotate(mol_idx=1, theta=180, axis='z')
        walton.rotate(mol_idx=1, theta=90, axis='y')
        walton()
        self.smiles_assertEqual('CC(C)(N)O', walton.merged)

    def test_tetraherdron_with_atom_change(self):
        walton = Walton.from_smiles(methylether='COC', methylether2='COC')
        for i in range(2):
            walton.flatten_trio(mol_idx=i, atom_idcs=(0, 1, 2))
            walton.atom_to_origin(mol_idx=i, atom_idx=1)
        walton.rotate(mol_idx=i, theta=180, axis='z')
        walton.rotate(mol_idx=i, theta=90, axis='y')
        walton()
        self.smiles_assertEqual('CC(C)(C)C', walton.merged)


if __name__ == '__main__':
    unittest.main()
