import os
from abc import abstractmethod
from itertools import chain
from typing import Union, Tuple, List, Dict

from examples.protocols_mergeCombineBase import Protocol_mergeCombineBase
from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.steps.hitsPreprocess_base import HitsPreprocess_base
from fragmenstein.protocols.steps.hitsPreprocess_fragmentationBrics import HitsPreprocess_fragmentationBRICS
from fragmenstein.protocols.steps.hitsPreprocess_permutations import HitsPreprocess_permutations
from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.utils.config_manager import ConfigManager


class Protocol_combineBase(Protocol_mergeCombineBase):

    def initialize(self, hit_ids,  preprocess_mode, template, template_xchemId, templates_dir=None,
                   template_pattern=None, max_attemps=None,
                   random_seed=None, *args, **kwargs):

        self._hit_ids = hit_ids
        self.preprocess_mode = preprocess_mode

        self.template_info_dict = dict(template=template,
                                       templates_dir=templates_dir,
                                       template_pattern=template_pattern)
        self.template_xchemId = template_xchemId

        self.max_attemps = max_attemps
        self.random_seed = random_seed if random_seed else ConfigManager.RANDOM_SEED_PERMUT


    @property
    def sdf_outname(self):
        return os.path.join(self.output_dir, ",".join(self.hit_ids) + ".sdf")

    @property
    @abstractmethod
    def permutations_size(self):
        raise NotImplementedError()

    @property
    @abstractmethod
    def combiner_class(self):
        raise NotImplementedError()

    @property
    def combinations_instead_permutations(self):
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

    def _permutate(self):

        if self.permutations_size:
            assert self.permutations_size >1 and self.permutations_size <= len(self.hit_ids), "Error, the permutations size should be between 2 and num_hits"

            proprocesser = HitsPreprocess_permutations(self.fragments, random_seed=self.random_seed)
            fragsCombin_iter = proprocesser.yield_combinations(max_num_elems=self.permutations_size, take_n_random=self.max_attemps,
                                                               combinations_instead_permutations = self.combinations_instead_permutations)
        else:
            fragsCombin_iter = [self.fragments]

        return fragsCombin_iter

    def _getFragmentatorClass(self):

        if self.preprocess_mode:
            print("preprocessing fragments")
            if self.preprocess_mode=="BRICS_decomposition":
                initialize_preprocesser = lambda  fragments: HitsPreprocess_fragmentationBRICS(fragments, random_seed=self.random_seed)
            elif self.preprocess_mode=="SMARTS_decomposition":
                raise NotImplementedError()
            else:
                raise ValueError("Error, not implemented option preprocess_mode=%s"%self.preprocess_mode)
        else:
            return None

        return initialize_preprocesser

    def preprocess_fragments(self) -> Tuple[Tuple[List[Compound], Dict[str, str]]]:

        fragsCombin_iter = self._permutate()
        initialize_preprocesser = self._getFragmentatorClass()

        bitId_to_molId  = {}
        if initialize_preprocesser:
            list_of_iters = []
            for fragGroup in fragsCombin_iter:
                preprocesser = initialize_preprocesser(fragGroup)
                list_of_iters.append(  preprocesser.yield_combinations(take_n_random=self.max_attemps) )
                bitId_to_molId.update( preprocesser.bitId_to_molId )
            fragsCombin_iter = chain.from_iterable( list_of_iters )

        fragsCombin_iter = list (fragsCombin_iter )
        return fragsCombin_iter, bitId_to_molId


    def checkTemplateArgs(self):

        if self.template_info_dict["template"] is None:
            assert self.template_xchemId is  None, "Error, template_xchemId cannot be provided if templates_dir provided "
            assert self.template_info_dict["template"] is  None, "Error, template is incompatible with templates_dir "
            assert self.template_info_dict["template_pattern"] is not None, "Error, template_pattern should be provided if templates_dir provided "
        else:

            assert self.template_xchemId is  not None, "Error, template_xchemId should be provided if template provided "
            assert self.template_info_dict["templates_dir"] is  None, "Error, template_pattern is incompatible with template "
            assert self.template_info_dict["template_pattern"] is  None, "Error, template_pattern is incompatible with template "

    def compute(self):
        print("computing")

        fragsCombin_iter, fragId_to_oriFragId = self.preprocess_fragments()
        self.checkTemplateArgs()

        combiner = self.combiner_class(output_path=self.wdir_fragmenstein, use_dask =ConfigManager.N_CPUS > 1,
                                       ** self.template_info_dict)

        results = combiner.applyCombine(fragsCombin_iter)

        for comp in results:
            if fragId_to_oriFragId:
                frags = comp.ref_molIds
                comp.ref_molIds = [fragId_to_oriFragId.get(frag, None) for frag in frags]

        return  results


    @classmethod
    def generateCmdParser(cls, prog, description):
        from fragmenstein.utils.cmd_parser import ArgumentParser
        parser = ArgumentParser(prog=prog, description=description)
        parser.add_argument("-i", "--data_root_dir", type=str,
                            help="The Xchem root dir for data, typically target_name/aligned/ ", required=True)
        parser.add_argument("-f", "--hit_ids", type=str, nargs="+", help="The hits ids to use  in the form x0020",
                            required=True)
        parser.add_argument("-o", "--output_dir", type=str, help="The directory where results will be saved",
                            required=True)

        parser.add_argument("-t", "--template", type=str,
                            help="The path to a template pdb. If not provided, the first hit found in --templates_dir"
                                 "would be used instead",
                            required=False)
        parser.add_argument("-x", "--template_xchemId", type=str,
                            help="The xchem id that would be used for reference pdb in Fragalysis if --template used", required=False)

        parser.add_argument("-d", "--templates_dir", type=str,
                            help="The directory where templates would be look for if not --template", required=False)

        parser.add_argument("--template_pattern", type=str,
                            help="The regex pattern of a template to find them in --templates_dir. "
                                 "Default: '%(default)s'", required=False, default= Xchem_info.unboundPdb_id_pattern)


        parser.add_argument("-s", "--permutations_size", type=int, default=None,
                            help="The number of fragments to process together. Default all. Min 2")

        parser.add_argument("-p", "--preprocess_mode", type=str,
                            choices=["BRICS_decomposition", "SMARTS_decomposition"], default=None,
                            help="How to preprocess the fragments. XXX_decomposition breaks the fragments into bits.",
                            required=False)

        parser.add_argument("-m", "--max_attemps", type=int, help="The number of maximun random attempts",
                            required=False, default=None)

        return parser
