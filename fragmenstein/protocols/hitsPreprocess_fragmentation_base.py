import random
from collections import defaultdict
from itertools import combinations, chain, product

from rdkit import Chem
from rdkit.Chem.Descriptors import HeavyAtomMolWt

from fragmenstein.external.smarts_fragmentation.reaction_fragmentation import ReactionFragmentation
from fragmenstein.protocols.hitsPreprocess_base import HitsPreprocess_base
from fragmenstein.utils.compound import Compound
from fragmenstein.utils.modify_molecules import change_atom_type


class HitsPreprocess_fragmentation_base(HitsPreprocess_base) :


    TEMPORAL_DUMMY_ATOM ="[Xe]"
    MAX_NUMBER_TO_COMBINE = 3

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

        self.yield_possible_molFrags = fragmentation_function

        self._all_fragmentations = None # ( primary_compound: [ fragOption1 ... fragOptionN] )
                                            # fragOption_i : [  frag1, frag2... ]
        self.bitId_to_molId = None

        super().__init__(original_fragments=original_fragments, random_seed=random_seed)

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
            bitId_to_molId[mol_id] = mol_id
            molId_to_bitLists_dict[mol_id].append([_mol])
            mol = change_atom_type(_mol, initial_symbol= '*', final_symbol=HitsPreprocess_fragmentation_base.TEMPORAL_DUMMY_ATOM)
            for i, split_option in  enumerate(self.yield_possible_molFrags(mol)):
                bits = sorted(split_option, key=lambda x: HeavyAtomMolWt(x))
                molId_to_bitLists_dict[mol_id].append([])
                for j, bit in enumerate(bits):
                    bitId = mol_id + "_%db%d" % (i,j)
                    bit = Chem.DeleteSubstructs(bit, Chem.MolFromSmiles('*'))

                    bit = change_atom_type(bit, initial_symbol=HitsPreprocess_fragmentation_base.TEMPORAL_DUMMY_ATOM, final_symbol='*')
                    bit = Compound(bit)
                    bit.molId = bitId
                    bit.parents = [_mol]
                    bitId_to_molId[bitId] = mol_id
                    molId_to_bitLists_dict[mol_id][-1].append( bit )

        return molId_to_bitLists_dict, bitId_to_molId


    def yield_combinations(self, min_num_fragments= 2, take_n_random=None):
        '''
        yields [compound1, compound2... compoundN]
        :param min_num_fragments:
        :param take_n_random:
        :return:
        '''
        assert len(self.original_fragments_dict) <= HitsPreprocess_fragmentation_base.MAX_NUMBER_TO_COMBINE, "Error, yield combinations scales with O(2**N_bits)**N_FRAGS, " \
                                                                   "so the number of fragments to consider has been limited to %d" % HitsPreprocess_fragmentation_base.MAX_NUMBER_TO_COMBINE

        full_compounds = list( self.original_fragments_dict.values() )

        try:
            final_fragments = self._yield_combinations(min_num_fragments)
            final_fragments = self.take_random_from_iterator(final_fragments, take_n_random)
        except ValueError:
            final_fragments = self._yield_combinations(min_num_fragments)

        final_fragments = chain.from_iterable([[full_compounds], final_fragments])
        return final_fragments

    def _yield_combinations(self, min_num_fragments=2):

        bitDecompositions_perCompound = self.broken_fragments.values() # [compound1_all_options, compound2_all_options, ...]

        oneDecomposition_perCompound_list =  list( product( *bitDecompositions_perCompound ))

        final_fragments = []
        for oneDecomposition_perCompound in oneDecomposition_perCompound_list:
            oneDecomposition_perCompound = list(chain.from_iterable(oneDecomposition_perCompound))
            enum_bits = self._powerset(oneDecomposition_perCompound, min_num_elements= min_num_fragments, include_full=False)
            enum_bits = filter( lambda bits_list: len(set( bit.primitiveId for bit in bits_list))>=min_num_fragments, enum_bits)
            final_fragments.append(enum_bits)

        final_fragments = chain.from_iterable( final_fragments)
        return  final_fragments

    def getOrinalFragmentId(self, bit):
        if isinstance(bit, Compound):
            bitId = bit.molId
        elif isinstance(bit, str):
            bitId = bit
        return self.bitId_to_molId[bitId]