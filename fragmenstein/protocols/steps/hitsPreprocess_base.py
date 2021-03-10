import random
from abc import abstractmethod
from itertools import chain, permutations
from itertools import combinations

from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.steps.adapt_input import InputAdapter


class HitsPreprocess_base(InputAdapter):


    @staticmethod
    def get_examples_init_params():

        input_mols = [("biphenyl", "c1ccccc1-c2ccccc2"),
                      ("ADP", "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N"),
                      ("PTC", "C1=CC=C(C=C1)NC(=S)N")]

        def create( name_smi):
            c = Compound.MolFromSmiles(name_smi[1])
            c.molId = name_smi[0]
            return c

        input_mols = list(map( create , input_mols))
        return dict( original_fragments = input_mols )

    def __init__(self,  original_fragments, random_seed=None):

        self.original_fragments_dict = self.adapt_dict_or_compoundsList( original_fragments )
        self.random_seed = random_seed


    def getOrinalFragmentId(self, mol):
        '''
        This should be overwriten by all class that either fragments or combines fragments
        :param mol:
        :return:
        '''

        if isinstance(mol, Compound):
            molId = mol.molId
        elif isinstance(mol, str):
            molId = mol
        return molId

    @classmethod
    def take_random_from_iterator(cls, iterable, take_n_random, random_seed=None):

        if take_n_random:
            if random_seed:
                random.seed(random_seed)
            iterable = cls.iter_sample_fast(iterable, take_n_random)
            if random_seed:
                random.seed(None)
            # iterable = list( iterable)
            return  iterable
        else:
            return iterable

    def _take_random_from_iterator(self, iterable, take_n_random):

        return type(self).take_random_from_iterator(iterable, take_n_random, self.random_seed)

    @abstractmethod
    def yield_combinations(self):
        '''
            yields [compound1, compound2... compoundN] derived from self.original_fragments_dict
        '''
        raise NotImplementedError()

    @classmethod
    def _powerset(cls, iterElems, min_num_elements=2, max_num_emements=None, include_full=False, combinations_instead_permutations=True):
        s = list(iterElems)
        if max_num_emements is None:
            last = len(s)+1 if include_full else len(s)
        else:
            last = max_num_emements+1
        func = combinations if combinations_instead_permutations else permutations
        return chain.from_iterable(func(s, r) for r in range(min_num_elements, last))



    @classmethod
    def iter_sample_fast(cls, iterable, samplesize):

        '''
        https://stackoverflow.com/questions/12581437/python-random-sample-with-a-generator-iterable-iterator
        '''
        results = []
        iterator = iter(iterable)
        # Fill in the first samplesize elements:
        try:
            for _ in range(samplesize):
                results.append(next(iterator))
        except StopIteration:
            raise ValueError("Sample larger than population.")
        random.shuffle(results)  # Randomize their positions
        for i, v in enumerate(iterator, samplesize):
            r = random.randint(0, i)
            if r < samplesize:
                results[r] = v  # at a decreasing rate, replace random items
        return results
