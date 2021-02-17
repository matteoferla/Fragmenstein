from typing import List

from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem.Lipinski import RotatableBondSmarts

def split_mol_to_brics_bits( mol: Chem.Mol, delete_dummy=False) -> List[Chem.Mol]:
    '''
    Split a molecule into a list of smaller molecules.

    :param mol. Chem.Mol object to be broken up into fragments by breaking rotable bonds
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
