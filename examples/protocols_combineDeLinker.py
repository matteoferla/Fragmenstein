from examples.protocols_combineBase import Protocol_combineBase
from fragmenstein.protocols.steps.combineMerge_DeLinkerDefault import CombineMerge_DeLinkerDefault

class Protocol_combineFragmenstein(Protocol_combineBase):

    @property
    def permutations_size(self):
        return 2

    @property
    def combinations_instead_permutations(self):
        return True

    @property
    def combiner_class(self):
        return CombineMerge_DeLinkerDefault


if __name__ == "__main__":

    parser = Protocol_combineFragmenstein.generateCmdParser(prog="protocols_combineDeLinker",
                                                            description="Combines molecules using fragmentation and DeLinker")


    args = vars(parser.parse_args())
    print(args)
    Protocol_combineFragmenstein.main(**args)
    print("\nmain DONE!\n")

    #TODO: modify fragalysis writer to accept method=DeLInker

    '''

N_CPUS=1 python -m examples.protocols_combineDeLinker -d ~/oxford/myProjects/diamondCovid/data/nsp13/aligned -f x0176_0B x0246_0B x0438_0B -o ~/oxford/tools/Fragmenstein/output -m 10 -p BRICS_decomposition


python -m fragmenstein.external.condor_queue.send_to_condor --env_vars EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json DASK_WORKER_MEMORY=4GB --ncpus 8 "/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python -m examples.protocols_combineDeLinker --n_cpus 8 -d /data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/aligned -o /data/xchem-fragalysis/sanchezg/oxford/tools/Fragmenstein/output_nsp13_Site7_mixed -f x0116_0B x0309_0B x4094_0B"

    '''