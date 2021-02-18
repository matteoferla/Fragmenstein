from .reactions_smarts import  reaction_smarts
from rdkit import Chem
from rdkit.Chem import AllChem

class ReactionFragmentation():

    def __init__(self, reactions_dict = reaction_smarts ):
        self.reactions_dict = reactions_dict

    def yield_reactant_decomposition(self, *mol_reactants):
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
                    for product in reaction:
                        # print(Chem.MolToSmiles(product))
                        products.append(product)
                    product_smiles= ",".join( sorted( map(Chem.MolToSmiles, products)))
                    if not product_smiles in all_decompositions:
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