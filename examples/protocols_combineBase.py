import os
from abc import abstractmethod
from typing import Tuple, List, Dict

from examples.protocols_mergeCombineBase import Protocol_mergeCombineBase
from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.utils.config_manager import ConfigManager


class Protocol_combineBase(Protocol_mergeCombineBase):

    def initialize(self, hit_ids,  preprocess_mode, data_root_dir, template, template_xchemId, templates_dir=None,
                   template_pattern=None, max_attemps=None, random_seed=None, protocol_verbose=False, *args, **kwargs):

        self._hit_ids = hit_ids
        self.data_root_dir = data_root_dir
        self.preprocess_mode = preprocess_mode

        self.template_info_dict = dict(template=template,
                                       templates_dir=templates_dir,
                                       template_pattern=template_pattern)
        self.template_xchemId = template_xchemId

        self.max_attemps = max_attemps
        self.random_seed = random_seed if random_seed else ConfigManager.RANDOM_SEED_PERMUT
        self.verbose = protocol_verbose

    @property
    def sdf_outname(self):
        return os.path.join(self.output_dir, ",".join(self.hit_ids) + ".sdf")


    @property
    @abstractmethod
    def combiner_class(self):
        raise NotImplementedError()


    @property
    def hit_ids(self):
        return self._hit_ids

    @hit_ids.setter
    def hit_ids(self, value):
        self._hit_ids = value

    @property
    def sdf_outname(self):
        return os.path.join(self.output_dir, ",".join(self.hit_ids) + ".sdf")


    @abstractmethod
    def preprocess_fragments(self) -> Tuple[Tuple[List[Compound], Dict[str, str]]]:
        raise NotImplementedError()


    def checkTemplateArgs(self):

        if self.template_info_dict["template"] is None: #No template provided. Then use --templates_dir or --data_root_dir
            assert self.template_xchemId is None, "Error, template_xchemId cannot be provided if templates_dir provided "

            if not self.template_info_dict["templates_dir"]:
                self.template_info_dict["templates_dir"] = self.data_root_dir

            assert self.template_info_dict["template"] is  None, "Error, template is incompatible with templates_dir "
            assert self.template_info_dict["template_pattern"] is not None, "Error, template_pattern should only be provided if templates_dir provided "

        else: #a template filename was provided

            assert self.template_xchemId is not None, "Error, template_xchemId should be provided if template provided "
            assert self.template_info_dict["templates_dir"] is  None, "Error, template_pattern is incompatible with template "
            assert self.template_info_dict["template_pattern"] is  None, "Error, template_pattern is incompatible with template "

    def compute(self):
        print("computing enumeration", flush=True)

        fragsCombin_iter, fragId_to_oriFragId = self.preprocess_fragments()
        self.checkTemplateArgs()

        combiner = self.combiner_class(output_path=self.wdir_enumeration, use_dask =ConfigManager.N_CPUS > 1,
                                       verbose = self.verbose, ** self.template_info_dict)

        results = combiner.applyCombine(fragsCombin_iter)

        for comp in results:
            if fragId_to_oriFragId:
                frags = comp.ref_molIds
                comp.ref_molIds = [fragId_to_oriFragId.get(frag, None) for frag in frags]

        return  results


    @classmethod
    def generateCmdParser(cls, prog, description):
        parser = Protocol_mergeCombineBase.generateCmdParser(prog, description)

        #TODO: skip_enummeration

        parser.add_argument("-i", "--data_root_dir", type=str,
                            help="The Xchem root dir for data, typically target_name/aligned/ ", required=True)
        parser.add_argument("-f", "--hit_ids", type=str, nargs="+", help="The hits ids to use  in the form x0020",
                            required=True)

        parser.add_argument("-t", "--template", type=str,
                            help="The path to a template pdb. If not provided, the first hit found in --templates_dir"
                                 "would be used instead. Default value same as --i",
                            required=False)
        parser.add_argument("-x", "--template_xchemId", type=str,
                            help="The xchem id that would be used for reference pdb in Fragalysis if --template used", required=False)

        parser.add_argument("-s", "--permutations_size", type=int, default=None,
                            help="The number of fragments to process together. Default all. Min 2")

        parser.add_argument("-p", "--preprocess_mode", type=str,
                            choices=["BRICS_decomposition", "SMARTS_decomposition"], default=None,
                            help="How to preprocess the fragments. XXX_decomposition breaks the fragments into bits.",
                            required=False)

        parser.add_argument("-m", "--max_attemps", type=int, help="The number of maximun random attempts",
                            required=False, default=None)

        return parser
