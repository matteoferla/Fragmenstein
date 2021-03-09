from rdkit import Chem

from fragmenstein.protocols.steps.hitsPreprocess_base import HitsPreprocess_base
from fragmenstein.utils.config_manager import ConfigManager


class HitsPreprocess_permutations(HitsPreprocess_base):


    def __init__(self, original_fragments, random_seed=None, *args, **kwargs):

        super().__init__(original_fragments=original_fragments, random_seed=random_seed)


    def yield_combinations(self,  max_num_elems=ConfigManager.COMBINE_PERMUTATIONS_MAX_NUM_ELEMENTS, take_n_random=None,
                           combinations_instead_permutations = True):

        get_results = lambda : self._powerset( self.original_fragments_dict.values() , min_num_elements=2,
                                               max_num_emements=max_num_elems,
                                               combinations_instead_permutations = combinations_instead_permutations)
        results = get_results()
        try:
            results = self._take_random_from_iterator(results, take_n_random)
        except ValueError:
            results = get_results()
        return results


def test_preprocessPermutations():

    fragmentator = HitsPreprocess_permutations( ** HitsPreprocess_permutations.get_examples_init_params()  )
    combs = fragmentator.yield_combinations()
    print(next( combs ) )
    for  option in combs:
        print( list( map(lambda x: (x.molId, x.primitiveId, Chem.MolToSmiles(x)), option) ) )


if __name__ =="__main__":
    test_preprocessPermutations()

    '''

python -m fragmenstein.protocols.hitsPreprocess_permutations

    '''