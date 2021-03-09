from itertools import chain
from typing import Tuple, List, Dict

from examples.protocols_combineBase import Protocol_combineBase
from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.steps.combineMerge_fragmensteinDefault import CombineMerge_FragmensteinDefault
from fragmenstein.protocols.steps.hitsPreprocess_base import HitsPreprocess_base
from fragmenstein.protocols.steps.hitsPreprocess_fragmentationBrics import HitsPreprocess_fragmentationBRICS
from fragmenstein.protocols.steps.hitsPreprocess_permutations import HitsPreprocess_permutations


class Protocol_combineFragmenstein(Protocol_combineBase):

    def initialize(self, permutations_size= None, *args, **kwargs): #added permutations_size to the builder

        super().initialize(*args, **kwargs)
        self.permutations_size = permutations_size


    @property
    def combiner_class(self):
        return CombineMerge_FragmensteinDefault

    def _permutate(self, fragments=None):


        if fragments == None:
            fragments = self.fragments

        if self.permutations_size:
            assert self.permutations_size >1 and self.permutations_size <= len(self.hit_ids), "Error, the permutations size should be between 2 and num_hits"

            proprocesser = HitsPreprocess_permutations(fragments, random_seed=self.random_seed)
            fragsCombin_iter = proprocesser.yield_combinations(max_num_elems=self.permutations_size, take_n_random=self.max_attemps,
                                                               combinations_instead_permutations = False)
        else:
            fragsCombin_iter = [fragments]

        return fragsCombin_iter

    def _getFragmentatorFactory(self, *args, **kwargs):

        initialize_preprocesser = None
        if self.preprocess_mode:
            print("preprocessing fragments")
            if self.preprocess_mode=="BRICS_decomposition":
                initialize_preprocesser = lambda  fragments: HitsPreprocess_fragmentationBRICS(fragments, *args, **kwargs,
                                                                                               random_seed=self.random_seed)
            elif self.preprocess_mode=="SMARTS_decomposition":
                raise NotImplementedError()
            else:
                raise ValueError("Error, not implemented option preprocess_mode=%s"%self.preprocess_mode)
        return initialize_preprocesser

    def preprocess_fragments(self) -> Tuple[Tuple[List[Compound], Dict[str, str]]]:


        fragsCombin_iter = self._permutate()
        initialize_preprocesser = self._getFragmentatorFactory()

        bitId_to_molId  = {}
        if initialize_preprocesser:
            list_of_iters = []
            for fragGroup in fragsCombin_iter:
                preprocesser = initialize_preprocesser(fragGroup)
                list_of_iters.append(  preprocesser.yield_combinations(take_n_random=self.max_attemps) )
                bitId_to_molId.update( preprocesser.bitId_to_molId )
            fragsCombin_iter = chain.from_iterable( list_of_iters )

        fragsCombin_iter = fragsCombin_iter

        fragsCombin_iter = HitsPreprocess_base.take_random_from_iterator( fragsCombin_iter, self.max_attemps)

        return fragsCombin_iter, bitId_to_molId


if __name__ == "__main__":

    parser = Protocol_combineFragmenstein.generateCmdParser(prog="protocol_combineFragmenstein", description="Combines molecules using Fragmenstein amd different fragmentations")

    # parser.add_argument( "--merging_mode", nargs=None, choices=["full", "partial", "none", "none_permissive", "off"], default="none_permissive",
    #                      help="See https://github.com/matteoferla/Fragmenstein/blob/master/documentation/monster/monster.md")

    args =vars( parser.parse_args())
    print(args)
    Protocol_combineFragmenstein.main( ** args)
    print("\nmain DONE!\n")

    '''

N_CPUS=1 python -m examples.protocols_combineFragmenstein -i ~/oxford/myProjects/diamondCovid/data/nsp13/aligned -f x0176_0B x0246_0B x0438_0B -o ~/oxford/tools/Fragmenstein/output -m 10 -p BRICS_decomposition


python -m fragmenstein.external.condor_queue.send_to_condor --env_vars EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json DASK_WORKER_MEMORY=4GB --ncpus 8 "/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python -m examples.protocols_combineFragmenstein --n_cpus 8 -d /data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/aligned -o /data/xchem-fragalysis/sanchezg/oxford/tools/Fragmenstein/output_nsp13_Site7_mixed -f x0116_0B x0309_0B x4094_0B"

    '''