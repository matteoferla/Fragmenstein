import random
from collections import defaultdict
from itertools import combinations, chain, product

from rdkit import Chem
from rdkit.Chem.Descriptors import HeavyAtomMolWt

from fragmenstein.external.smarts_fragmentation.reaction_fragmentation import ReactionFragmentation
from fragmenstein.protocols.adapt_input import InputAddapter
from fragmenstein.utils.compound import Compound
from fragmenstein.utils.modify_molecules import change_atom_type
from fragmenstein.utils.sample_iterator import iter_sample_fast


class Fragmentator(InputAddapter):


    TEMPORAL_DUMMY_ATOM ="[Xe]"
    MAX_NUMBER_TO_COMBINE = 4

    @staticmethod
    def get_examples_init_params():

        input_mols = [("biphenyl", "c1ccccc1-c2ccccc2"),
                      ("ADP", "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N")]

        def create( name_smi):
            c = Compound.MolFromSmiles(name_smi[1])
            c.molId = name_smi[0]
            return c

        input_mols = list(map( create , input_mols))
        return dict( fragments_dict = input_mols )

    def __init__(self, original_fragments, fragmentation_function, min_num_fragments= 2, take_n_random=None, random_seed=None,
                 *args, **kwargs):


        # if self.fragmentation_mode == Fragmentator.SMARTS_MODE:
        #     reactFragmentator = ReactionFragmentation(*args, **kwargs)
        #     self.yield_possible_molFrags = lambda mol : reactFragmentator.yield_reactant_decomposition(mol)
        # elif self.fragmentation_mode == Fragmentator.BRICS_MODE:
        #     self.yield_possible_molFrags = lambda mol: (split_mol_to_brics_bits(mol),)
        # else:
        #     raise ValueError("Error, not valid fragmentation_mode (%s) "%(fragmentation_mode, ) )


        self.min_num_fragments = min_num_fragments
        self.take_n_random = take_n_random

        self.random_seed = random_seed

        self.original_fragments_dict = self.adapt_dict_or_compoundsList( original_fragments )

        self.yield_possible_molFrags = fragmentation_function

        self._all_fragmentations = None # ( primary_compound: [ fragOption1 ... fragOptionN] )
                                            # fragOption_i : [  frag1, frag2... ]
        self.bitId_to_molId = None
    @property
    def broken_fragments(self):
        if self._all_fragmentations is None:
            self._all_fragmentations, self.bitId_to_molId = self.fragment_molsList( self.original_fragments_dict )
        return self._all_fragmentations




    def fragment_molsList(self, input_mols):
        '''

        :param input_mols: Molecules to fragment.2 Options
                            1) A Dict like object molId -> Compound
                            2) A list of Compounds with molId set
        :return:
        '''
        try:
            input_mols =  input_mols.items()
        except AttributeError:
            input_mols = [ (comp.molId,  comp) for comp in input_mols ]

        molId_to_bitLists_dict = defaultdict(list)
        bitId_to_molId = {}
        for mol_id, _mol in input_mols:
            mol = change_atom_type(_mol, initial_symbol= '*', final_symbol=Fragmentator.TEMPORAL_DUMMY_ATOM)
            for i, split_option in  enumerate(self.yield_possible_molFrags(mol)):
                bits = sorted(split_option, key=lambda x: HeavyAtomMolWt(x))
                molId_to_bitLists_dict[mol_id].append([])
                for j, bit in enumerate(bits):
                    bitId = mol_id + "_%d-%d" % (i,j)
                    bit = Chem.DeleteSubstructs(bit, Chem.MolFromSmiles('*'))

                    bit = change_atom_type(bit, initial_symbol=Fragmentator.TEMPORAL_DUMMY_ATOM, final_symbol= '*' )
                    bit = Compound(bit)
                    bit.molId = bitId
                    bit.parents = [_mol]
                    bitId_to_molId[bitId] = mol_id
                    molId_to_bitLists_dict[mol_id][-1].append( bit )

        return molId_to_bitLists_dict, bitId_to_molId


    def yield_bits_combinations(self, min_num_fragments= 2, take_n_random=None ):
        assert  len(self.original_fragments_dict) < Fragmentator.MAX_NUMBER_TO_COMBINE, "Error, yield combinations scales with O(2**N_bits)**N_FRAGS, " \
                                                                   "so the number of fragments to consider has been limited to %d"%Fragmentator.MAX_NUMBER_TO_COMBINE

        bitDecompositions_perCompound = self.broken_fragments.values() # [compound1_all_options, compound2_all_options, ...]

        oneDecomposition_perCompound_list =  list(product( *bitDecompositions_perCompound ))

        final_fragments = []
        full_compounds = list( self.original_fragments_dict.values() )

        for oneDecomposition_perCompound in oneDecomposition_perCompound_list:
            oneDecomposition_perCompound = list(chain.from_iterable(oneDecomposition_perCompound))
            enum_bits = self._powerset(oneDecomposition_perCompound, min_num_fragments= min_num_fragments)
            enum_bits = filter( lambda bits_list: len(set( bit.primitiveId for bit in bits_list))>=min_num_fragments, enum_bits)
            final_fragments.append(enum_bits)

        final_fragments = chain.from_iterable( final_fragments )

        if take_n_random:
            if self.random_seed:
                random.seed( self.random_seed )
            final_fragments = iter_sample_fast(final_fragments, take_n_random)
            if self.random_seed:
                random.seed(None)

        final_fragments = chain.from_iterable([ [full_compounds], final_fragments ])

        return final_fragments

    def _powerset(self, iterElems, min_num_fragments=2, include_full=False):
        s = list(iterElems)
        last = len(s)+1 if include_full else len(s)
        return chain.from_iterable(combinations(s, r) for r in range(min_num_fragments, last))

    def getOrinalFragmentId(self, bit):
        if isinstance(bit, Compound):
            bitId = bit.molId
        elif isinstance(bit, str):
            bitId = bit
        return self.bitId_to_molId[bitId]