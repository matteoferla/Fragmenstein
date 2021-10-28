from rdkit import Chem

from fragmenstein.external.smarts_fragmentation.reaction_fragmentation import ReactionFragmentation
from fragmenstein.pipelines.protocols.steps.hitsPreprocess_fragmentationBase import HitsPreprocess_fragmentationBase





class HitsPreprocess_fragmentationSMARTS(HitsPreprocess_fragmentationBase):


    def __init__(self, original_fragments, *args, **kwargs):

        fragmentator = ReactionFragmentation()
        super().__init__(original_fragments,  lambda mol: list(fragmentator.yield_reactant_decomposition(mol)), *args, **kwargs)

        raise NotImplementedError("Error")
        #TODO: find non matching atoms and replace with dummy



def test_fragmentBrics():
    raise NotImplementedError()
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