import os


from examples.protocols_FragmensteinBase import Protocol_FragmensteinBase
from fragmenstein.protocols.combineMerge_fragmensteinDefault import CombineMerge_FragmensteinDefault
from fragmenstein.protocols.hitsPreprocess_fragmentationBrics import HitsPreprocess_fragmentationBRICS
from fragmenstein.protocols.hitsPreprocess_permutations import HitsPreprocess_permutations

from fragmenstein.utils.config_manager import ConfigManager

class Protocol_combineFragmenstein(Protocol_FragmensteinBase):

    RANDOM_SEED =121

    def initialize(self, hit_ids,  preprocess_mode, template, template_xchemId, max_attemps, random_seed=None, *args, **kwargs):
        self.hit_ids = hit_ids
        self.preprocess_mode = preprocess_mode
        self.template = template
        self.template_xchemId = template_xchemId
        self.max_attemps = max_attemps
        self.random_seed = random_seed if random_seed else Protocol_combineFragmenstein.RANDOM_SEED

    @property
    def sdf_outname(self):
        return os.path.join(self.output_dir, ",".join(self.hit_ids) + ".sdf")

    def compute(self):
        print("computing")
        proprocesser = None
        if self.preprocess_mode:
            print("preprocessing fragments")
            if self.preprocess_mode=="BRICS_decomposition":
                proprocesser = HitsPreprocess_fragmentationBRICS(self.fragments, random_seed=self.random_seed )
                fragsCombin_iter =  proprocesser.yield_combinations(take_n_random=self.max_attemps)
            elif self.preprocess_mode=="SMARTS_decomposition":
                raise NotImplementedError()
            elif self.preprocess_mode == "permutations":
                proprocesser = HitsPreprocess_permutations(self.fragments, random_seed=self.random_seed )
                fragsCombin_iter =  proprocesser.yield_combinations(max_num_elems=3, take_n_random=self.max_attemps)
            else:
                raise ValueError("Error, not implemented option preprocess_mode=%s"%self.preprocess_mode)
        else:
            fragsCombin_iter = [self.fragments]

        if self.template is None:
            id_template_iter = self.data_loader.find_templates( filter_ids=[self.hit_ids[0]])
            template_id, template = list( id_template_iter)[0]
            if not self.template_xchemId:
                self.template_xchemId = template_id
        else:
            template = os.path.expanduser(self.template)
            assert self.template_xchemId is not None, "Error, template_xchemId should not be None if template is not None "

        combiner = CombineMerge_FragmensteinDefault(output_path=self.wdir_fragmenstein, template=template,
                                                    use_dask =ConfigManager.N_CPUS > 1)

        results = combiner.applyCombine(fragsCombin_iter)

        for comp in results:
            if proprocesser:
                frags = comp.ref_molIds
                comp.ref_molIds = [proprocesser.getOrinalFragmentId(frag) for frag in frags]

        return  results


if __name__ == "__main__":
    from fragmenstein.utils.cmd_parser import ArgumentParser
    parser = ArgumentParser(prog="protocol_bricsFragmenstein", description="Combines molecules using Fragmenstein amd different fragmentations")
    parser.add_argument("-d", "--data_root_dir", type=str, help="The Xchem root dir for data, typically target_name/aligned/ ", required=True)
    parser.add_argument("-f", "--hit_ids", type=str, nargs="+", help="The hits ids to use  in the form x0020", required=True)
    parser.add_argument("-o", "--output_dir", type=str, help="The directory where results will be saved", required=True)
    parser.add_argument("-t", "--template", type=str, help="The path to a template pdb. If not provided, the first hit would be used", required=False)
    parser.add_argument("-x", "--template_xchemId", type=str, help="The xchem id that would be used for reference pdb in fragalysis", required=False)
    parser.add_argument("-p", "--preprocess_mode", type=str, choices=["BRICS_decomposition", "permutations", "SMARTS_decomposition"], default= None,
                        help="How to preprocess the fragments. XXX_decomposition breaks the fragments into bits.", required=False)
    # parser.add_argument( "--merging_mode", nargs=None, choices=["full", "partial", "none", "none_permissive", "off"], default="none_permissive",
    #                      help="See https://github.com/matteoferla/Fragmenstein/blob/master/documentation/monster/monster.md")
    parser.add_argument("-m", "--max_attemps", type=int, help="The number of maximun random attempts", required=False, default=None)

    args =vars( parser.parse_args())
    print(args)
    Protocol_combineFragmenstein.main( ** args)
    print("\nmain DONE!\n")

    '''

N_CPUS=1 python -m examples.protocols_combineFragmenstein -d ~/oxford/myProjects/diamondCovid/data/nsp13/aligned -f x0176_0B x0246_0B x0438_0B -o ~/oxford/tools/Fragmenstein/output -m 10 -p BRICS_decomposition


python -m fragmenstein.external.condor_queue.send_to_condor --env_vars EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json DASK_WORKER_MEMORY=4GB --ncpus 8 "/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python -m examples.protocols_combineFragmenstein --n_cpus 8 -d /data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/aligned -o /data/xchem-fragalysis/sanchezg/oxford/tools/Fragmenstein/output_nsp13_Site7_mixed -f x0116_0B x0309_0B x4094_0B"

    '''