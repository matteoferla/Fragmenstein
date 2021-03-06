import os
import sys

from examples.protocols_mergeCombineBase import Protocol_mergeCombineBase
from fragmenstein.protocols.steps.combineMerge_fragmensteinDefault import CombineMerge_FragmensteinDefault
from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.utils.config_manager import ConfigManager


class Protocol_placeFragmenstein(Protocol_mergeCombineBase):


    def initialize(self, smile_fragIds_list, *args, **kwargs):
        self.smile_fragIds_list = smile_fragIds_list

    @property
    def sdf_outname(self):
        return os.path.join(self.output_dir, "placed_smiles.sdf")

    def compute(self):
        print("computing")
        fragmensteiner = CombineMerge_FragmensteinDefault(output_path=self.wdir_fragmenstein, template=None,
                                                          templates_dir=self.data_root_dir,
                                                          template_pattern=Xchem_info.unboundPdb_id_pattern,
                                                          merging_mode=self.merging_mode,
                                                          use_dask=ConfigManager.N_CPUS > 1)
        def getFragments(fragIds):
            return [self.data_loader.fragments_dict[fragId] for fragId in fragIds]

        results = fragmensteiner.applyPlace(((smi, getFragments(fragIds)) for smi, fragIds in smile_fragIds_list))
        return results

if __name__ == "__main__":
    from fragmenstein.utils.cmd_parser import ArgumentParser, argparse

    parser = ArgumentParser(prog="protocol_placeFragmenstein", description="places a smiles given the inspirational fragments using Fragmenstein")
    parser.add_argument("-d", "--data_root_dir", type=str, help="The Xchem root dir for data, typically target_name/aligned/ ", required=True)
    parser.add_argument("-i", "--input", nargs=None, type=argparse.FileType('r'), default=sys.stdin, help="Tab separated file with two columns: smiles and fragment_ids."
                                                                   " Fragment_ids are commma separated. E.g.:\n"
                                                                  "CCO  x0020_0B,x0029_0A", required=True)
    parser.add_argument( "--merging_mode", nargs=None, choices=["full", "partial", "none", "none_permissive", "off"], default="none_permissive",
                         help="See https://github.com/matteoferla/Fragmenstein/blob/master/documentation/monster/monster.md")
    parser.add_argument("-o", "--output_dir", type=str, help="The directory where results will be saved", required=True)

    args =vars( parser.parse_args())
    print(args)
    smile_fragIds_list = [ line.split() for line in args["input"].read().splitlines() ]
    smile_fragIds_list = [ (smi, fragIds.split(",")) for smi, fragIds in smile_fragIds_list ]
    Protocol_placeFragmenstein.main( smile_fragIds_list = smile_fragIds_list, ** args)
    print("\nmain DONE!\n")

    '''

echo -e "Cc1ccncc1NC(=O)C(C)c1cccc(Cl)c1  x0020_0B,x0029_0A
CC(=O)N1CCN(Cc2cccs2)CC1 x0020_0B" | N_CPUS=1 python -m examples.protocols_placeFragmenstein -i - -d ~/oxford/myProjects/diamondCovid/data/nsp13/aligned -o ~/oxford/tools/Fragmenstein/output


    '''