from typing import List

from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem.Lipinski import RotatableBondSmarts

from fragmenstein.protocols.steps.hitsPreprocess_fragmentationBase import HitsPreprocess_fragmentationBase


def split_mol_to_brics_bits(mol: Chem.Mol, delete_dummy=False) -> List[Chem.Mol]:
    '''
    Split a molecule into a list of smaller molecules.

    :param mol. Chem.Mol object to be broken up into fragments by breaking rotable bonds
    :param delete_dummy. Delete dummy atoms?
    :return:  a list of Chem.Mol objects that represent the bits in which input mol was broken.
    '''
    # find the rotatable bonds
    bonds = mol.GetSubstructMatches(RotatableBondSmarts)

    bonds = [((x, y), (0, 0)) for x, y in bonds]
    p = BRICS.BreakBRICSBonds(mol, bonds=bonds)
    mols = [mol for mol in Chem.GetMolFrags(p, asMols=True)]
    if delete_dummy:
        dummy_mol = Chem.MolFromSmarts('[#0]')
        mols = [Chem.DeleteSubstructs(mol, dummy_mol) for mol in mols]
    return mols


class HitsPreprocess_fragmentationBRICS(HitsPreprocess_fragmentationBase):


    def __init__(self, original_fragments, break_mode="default", *args, **kwargs):

        if break_mode == "default":
            frag_fun = self.get_all_atoms_frags
        elif break_mode =="binary":
            frag_fun = self.yield_binary_splits
        else:
            raise NotImplementedError("Error, no more break_modes implemented yet")
        super().__init__(original_fragments,  frag_fun, *args, **kwargs)


    def get_all_atoms_frags(self, mol):
        return (split_mol_to_brics_bits(mol),)

    def yield_binary_splits(self, mol):

        bonds = mol.GetSubstructMatches(RotatableBondSmarts)
        bonds = [((x, y), (0, 0)) for x, y in bonds]
        for bond in bonds:
            p = BRICS.BreakBRICSBonds(mol, bonds=[bond])
            mols = [mol for mol in Chem.GetMolFrags(p, asMols=True)]
            yield  mols



def test_fragmentBrics():
    print( [ Chem.MolToSmiles(mol) for mol in HitsPreprocess_fragmentationBRICS.get_examples_init_params()['original_fragments'] ])
    fragmentator = HitsPreprocess_fragmentationBRICS(** HitsPreprocess_fragmentationBRICS.get_examples_init_params())
    combs = fragmentator.yield_combinations()
    next( combs )
    for molId, fragmentation_options in fragmentator._all_fragmentations.items():
        print("--------- mol id:", molId, "---------" )
        for frag_option in fragmentation_options:
            print( list( map(lambda x: (x.molId, x.primitiveId, Chem.MolToSmiles(x)), frag_option) ) )

def test_fragmentBrics2():
    print( [ Chem.MolToSmiles(mol) for mol in HitsPreprocess_fragmentationBRICS.get_examples_init_params()['original_fragments'] ])
    fragmentator = HitsPreprocess_fragmentationBRICS(** HitsPreprocess_fragmentationBRICS.get_examples_init_params(), break_mode="binary")
    combs = fragmentator.yield_combinations()
    next( combs )
    for molId, fragmentation_options in fragmentator._all_fragmentations.items():
        print("--------- mol id:", molId, "---------" )
        for frag_option in fragmentation_options:
            print( list( map(lambda x: (x.molId, x.primitiveId, Chem.MolToSmiles(x)), frag_option) ) )

# test_fragmentBrics()

if __name__ =="__main__":
    test_fragmentBrics2()

    '''

python -m fragmenstein.protocols.steps.hitsPreprocess_fragmentationBrics

    '''