
import os
from typing import List

from rdkit import Chem

from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.steps.combineMerge_abstract import ErrorInComputation, CombineMerge_Base
from fragmenstein.protocols.steps.combineMerge_fragmensteinDefault import CombineMerge_FragmensteinDefault
from fragmenstein.protocols.steps.hitsPreprocess_fragmentationBrics import HitsPreprocess_fragmentationBRICS
from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.utils.config_manager import ConfigManager


class CombineMerge_DeLinkerDefault( CombineMerge_FragmensteinDefault  ):



    # @staticmethod
    # def get_examples_combine_params():
    #     data_dir = os.path.abspath(os.path.join(Xchem_info.examples_dir, "hit_mols"))
    #     fnames = [os.path.join(data_dir, "Mpro-x0678.mol"), os.path.join(data_dir, "Mpro-x0434.mol")]
    #     list_of_fragments = [Compound.MolFromFile(fname, fname.split("-")[-1].split(".")[0]) for fname in fnames]
    #     list_of_fragments = list(HitsPreprocess_fragmentationBRICS(list_of_fragments, random_seed=121).yield_combinations())
    #     print( [ [Chem.MolToSmiles(mol) for mol in enum ] for enum in list_of_fragments ])
    #     input( len(list_of_fragments) )
    #     return dict( list_of_fragments = list_of_fragments )

    def __init__(self, random_seed=None, gpu_id=None, number_of_generation_per_valid=10, *args, **kwargs):

        super().__init__( args, **kwargs)

        self.random_seed = random_seed

        from fragmenstein.external.DeLinker.DeLinkerWrapper import DeLinkerWrapper
        self.delinker = DeLinkerWrapper( number_of_generation_per_valid=number_of_generation_per_valid,
                         n_atomPairs_attemps=3, n_cores=1, gpu_id=gpu_id,
                         interactive=False, random_seed= self.random_seed)


    def tryOneGeneric(self, merge_id, templateFname, fragments: List[Compound], wdir, *args, **kwargs):

        assert len(fragments) == 2, "Error, DeLinker only works for pairs of compounds"
        #TODO: not working
        raise  NotImplementedError()
        proposed = self.delinker.link_molecule_pair(*fragments)

        placed_result = []
        for i, proposal in enumerate(proposed):
            smi = Chem.MolToSmiles( proposal )
            placed_result.append(
                           super().tryOneGeneric( merge_id+"_"+str(i), templateFname, fragments, wdir, smi)
            )
        return placed_result

    def _select_pieces(self, frags_option):

        print( frags_option)
        pass


def test_applyCombine():
    combiner = CombineMerge_DeLinkerDefault( **CombineMerge_DeLinkerDefault.get_examples_init_params(), use_dask=False)
    results = combiner.applyCombine( **CombineMerge_DeLinkerDefault.get_examples_combine_params())
    print("RESULTS applyCombine:")
    print( results)



if __name__ == "__main__":

    print("trying combine")
    test_applyCombine()

    '''

python -m fragmenstein.protocols.steps.combineMerge_DeLinkerDefault

    '''