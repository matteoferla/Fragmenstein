from examples.protocols_combineBase import Protocol_combineBase
from fragmenstein.protocols.steps.combineMerge_fragmensteinDefault import CombineMerge_FragmensteinDefault


class Protocol_combineFragmenstein(Protocol_combineBase):

    def initialize(self, permutations_size= None, *args, **kwargs): #added permutations_size to the builder

        super().initialize(*args, **kwargs)
        self._permutations_size = permutations_size

    @property
    def permutations_size(self):
        return self._permutations_size

    @property
    def combinations_instead_permutations(self):
        return False

    @property
    def combiner_class(self):
        return CombineMerge_FragmensteinDefault


if __name__ == "__main__":

    parser = Protocol_combineFragmenstein.generateCmdParser(prog="protocol_combineFragmenstein", description="Combines molecules using Fragmenstein amd different fragmentations")

    # parser.add_argument( "--merging_mode", nargs=None, choices=["full", "partial", "none", "none_permissive", "off"], default="none_permissive",
    #                      help="See https://github.com/matteoferla/Fragmenstein/blob/master/documentation/monster/monster.md")

    args =vars( parser.parse_args())
    print(args)
    Protocol_combineFragmenstein.main( ** args)
    print("\nmain DONE!\n")

    '''

N_CPUS=1 python -m examples.protocols_combineFragmenstein -d ~/oxford/myProjects/diamondCovid/data/nsp13/aligned -f x0176_0B x0246_0B x0438_0B -o ~/oxford/tools/Fragmenstein/output -m 10 -p BRICS_decomposition


python -m fragmenstein.external.condor_queue.send_to_condor --env_vars EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json DASK_WORKER_MEMORY=4GB --ncpus 8 "/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python -m examples.protocols_combineFragmenstein --n_cpus 8 -d /data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/aligned -o /data/xchem-fragalysis/sanchezg/oxford/tools/Fragmenstein/output_nsp13_Site7_mixed -f x0116_0B x0309_0B x4094_0B"

    '''