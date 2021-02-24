from typing import List

from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem.Lipinski import RotatableBondSmarts

from fragmenstein.protocols.fragmentation_base import Fragmentator


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


class Fragmentator_BRICS(Fragmentator):


    def __init__(self, fragments_dict, *args, **kwargs):
        super().__init__(fragments_dict,  lambda mol: (split_mol_to_brics_bits(mol),), *args, **kwargs)



def test_fragmentBrics():

    fragmentator = Fragmentator_BRICS( ** Fragmentator_BRICS.get_examples_init_params()  )
    combs = fragmentator.yield_bits_combinations()
    next( combs )
    for molId, fragmentation_options in fragmentator._all_fragmentations.items():
        print("--------- mol id:", molId, "---------" )
        for frag_option in fragmentation_options:
            print( list( map(lambda x: (x.molId, x.primitiveId, Chem.MolToSmiles(x)), frag_option) ) )


# test_fragmentBrics()

if __name__ =="__main__":
    test_fragmentBrics()

    '''

python -m fragmenstein.protocols.fragmentation_brics

    '''