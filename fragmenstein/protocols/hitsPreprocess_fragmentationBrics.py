from typing import List

from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem.Lipinski import RotatableBondSmarts

from fragmenstein.protocols.hitsPreprocess_fragmentation_base import HitsPreprocess_base


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


class HitsPreprocess_fragmentationBRICS(HitsPreprocess_base):


    def __init__(self, original_fragments, *args, **kwargs):
        super().__init__(original_fragments,  lambda mol: (split_mol_to_brics_bits(mol),), *args, **kwargs)



def test_fragmentBrics():

    fragmentator = HitsPreprocess_fragmentationBRICS(** HitsPreprocess_fragmentationBRICS.get_examples_init_params())
    combs = fragmentator.yield_combinations()
    next( combs )
    for molId, fragmentation_options in fragmentator._all_fragmentations.items():
        print("--------- mol id:", molId, "---------" )
        for frag_option in fragmentation_options:
            print( list( map(lambda x: (x.molId, x.primitiveId, Chem.MolToSmiles(x)), frag_option) ) )


# test_fragmentBrics()

if __name__ =="__main__":
    test_fragmentBrics()

    '''

python -m fragmenstein.protocols.hitsPreprocess_fragmentationBrics

    '''