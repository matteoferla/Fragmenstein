import itertools
from functools import lru_cache

from rdkit.Chem import rdMolHash
from rdkit.Chem.Descriptors import HeavyAtomMolWt

from .reactions_smarts import  REACTION_SMARTS
from rdkit import Chem
from rdkit.Chem import AllChem

class ReactionFragmentation():

    def __init__(self, reactions_dict = REACTION_SMARTS, min_MW_bit = 40):
        self.reactions_dict = reactions_dict
        self.min_MW_bit = min_MW_bit

    def yield_decompositions_one_option(self, *mol_reactants):
        # print( Chem.MolToSmiles(*mol_reactants))
        all_decompositions = set([])
        for reactType, reactions in self.reactions_dict.items():
            for reactSmart in  reactions:
                # print( reactSmart )
                try:
                    rxn = AllChem.ReactionFromSmarts(reactSmart)
                    ps = rxn.RunReactants(mol_reactants)
                except RuntimeError:
                    continue
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

    def _hashDecomp(self, decomp): #TODO: define sorting criterium

        # mol_hash = ""
        # for comp in sorted(decomp, key=lambda x: Chem.MolToSmiles(x)):
        #     mol_hash += rdMolHash.MolHash(comp, rdMolHash.HashFunction.names['AnonymousGraph']) + \
        #                 rdMolHash.MolHash(comp, rdMolHash.HashFunction.names['CanonicalSmiles'])

        smiles = map(Chem.MolToSmiles, decomp)
        mol_hash = "-".join( smiles )

        return mol_hash

    def yield_reactant_decomposition(self, mol):

        seen_solutions = set([])
        for decomp in self._get_reactant_decomposition(mol):
            mol_hash = self._hashDecomp(decomp)
            if mol_hash not in seen_solutions:
                seen_solutions.add( mol_hash)
                yield decomp


    @lru_cache(maxsize=None)
    def _get_reactant_decomposition(self, mol):
        '''

        :param mol_reactants: Chem.Mol objects
        :return:
        '''
        # print ( mol_reactants )


        decompositions = list(self.yield_decompositions_one_option( mol ))

        if len(decompositions)==0:
            return []
        else:
            final_decomps = [ (mol,) ]
            for ondeDecomp in decompositions:
                # print(">>>>>", ondeDecomp)
                left_decomposition, right_decomposition = ondeDecomp
                left_decomposition_second = list(self._get_reactant_decomposition(left_decomposition)) + [[left_decomposition]]
                right_decomposition_second = list(self._get_reactant_decomposition(right_decomposition)) + [[right_decomposition]]

                # print("+++++++", left_decomposition_second)
                # print("*******", right_decomposition_second)
                # print("-------", list(map(list, itertools.product(left_decomposition_second, right_decomposition_second))))

                final_decomps.extend(tuple(map(lambda prod: tuple(itertools.chain.from_iterable(prod)), itertools.product(left_decomposition_second, right_decomposition_second))))
            # print(final_decomps)
            return final_decomps

def reaction_fragmentation_test1():
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

def reaction_fragmentation_test2():
    query = "aaaaaaaaa0bbbbbbbbb0cccccc0dddd0ee"

    def yield_decompositions_one_option_dummy(self, query):
        indices = [i for i, x in enumerate(query) if x == "0"]
        for i in indices:
            yield list((query[:i], query[i+1:]))

    def _hashDecomp(self, decomp):
        return "-".join( sorted(decomp))

    ReactionFragmentation.yield_decompositions_one_option = yield_decompositions_one_option_dummy
    ReactionFragmentation._hashDecomp = _hashDecomp

    results = list(ReactionFragmentation().yield_reactant_decomposition(query))

    print(sorted(results))
    print("--------------")
    print(sorted(set(results)))

def reaction_fragmentation_test3():
    query = "aaaaaaaaa0bbbbbbbbb*cccccc*dddd0ee"

    def yield_decompositions_one_option_dummy(self, query):
        for character in ["0", "*"]:
            indices = [i for i, x in enumerate(query) if x == character]
            for i in indices:
                yield list((query[:i], query[i+1:]))

    def _hashDecomp(self, decomp):
        return "-".join( sorted(decomp))

    ReactionFragmentation.yield_decompositions_one_option = yield_decompositions_one_option_dummy
    ReactionFragmentation._hashDecomp = _hashDecomp

    results = list(ReactionFragmentation().yield_reactant_decomposition(query))

    print(sorted(results))
    print("--------------")
    print(sorted(set(results)))

def reaction_fragmentation_test4():

    query_smiles = "Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12" # 'CC(C(C1=CC(=C(C=C1)O)O)O)N(C)C(=O)OCC2=CC=CC=C2'
    query_mol = Chem.MolFromSmiles(query_smiles)
    from rdkit.Chem.rdDistGeom import EmbedMolecule
    EmbedMolecule(query_mol)
    results = list(ReactionFragmentation().yield_reactant_decomposition(query_mol))

    print(results)
    from matplotlib import pyplot as plt
    from rdkit.Chem import Draw
    for reacts in results:
        print(  [Chem.MolToSmiles(elem) for elem in reacts] )
        conf = reacts[0].GetConformer(0)
        print(conf)
        try:
            plt.imshow(Draw.MolsToGridImage( [query_mol] + list(reacts), molsPerRow=1))
            plt.show()
        except RuntimeError:
            continue

def main():
    pass
if __name__ == "__main__":
    main()


    '''

python -m fragmenstein.external.smarts_fragmentation.reaction_fragmentation

    '''