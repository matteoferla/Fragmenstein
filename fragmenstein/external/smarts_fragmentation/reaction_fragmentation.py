from rdkit.Chem.Descriptors import HeavyAtomMolWt

from .reactions_smarts import  REACTION_SMARTS
from rdkit import Chem
from rdkit.Chem import AllChem

class ReactionFragmentation():

    def __init__(self, reactions_dict = REACTION_SMARTS, min_MW_bit = 40):
        self.reactions_dict = reactions_dict
        self.min_MW_bit = min_MW_bit

    def yield_reactant_decomposition(self, *mol_reactants):
        '''

        :param mol_reactants: Chem.Mol objects
        :return:
        '''
        # print ( mol_reactants )

        all_decompositions = set([])
        for reactType, reactions in self.reactions_dict.items():
            for reactSmart in  reactions:
                # print( reactSmart )
                rxn = AllChem.ReactionFromSmarts(reactSmart)
                ps = rxn.RunReactants(mol_reactants)
                for i, reaction in enumerate(ps):
                    products = []
                    # print("reaction %d"%i)
                    ignore = False
                    for product in reaction:
                        # print(Chem.MolToSmiles(product))
                        if HeavyAtomMolWt(product) < self.min_MW_bit:
                            ignore=True
                            break
                        products.append(product)
                    product_smiles= ",".join( sorted( map(Chem.MolToSmiles, products)))
                    if not product_smiles in all_decompositions and not ignore:
                        all_decompositions.add(product_smiles)
                        yield products

def example():
    query_smiles = "COCC(=O)Nc1[n]c2ccccc2[s]1"
    query_mol = Chem.MolFromSmiles(query_smiles)
    results = ReactionFragmentation().yield_reactant_decomposition(query_mol)
    # print( [ [Chem.MolToSmiles(elem) for elem in reacts] for reacts in results] )

    from matplotlib import pyplot as plt
    from rdkit.Chem import Draw
    for reacts in results:
        print(  [Chem.MolToSmiles(elem) for elem in reacts] )
        plt.imshow(Draw.MolsToGridImage( [query_mol] + list(reacts), molsPerRow=1))
        plt.show()

if __name__ == "__main__":
    example()

    '''

python -m fragmenstein.external.smarts_fragmentation.reaction_fragmentation

    '''