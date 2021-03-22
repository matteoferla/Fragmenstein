import os
import pickle
import shutil
import tempfile
from abc import ABC, abstractmethod
from itertools import chain
from typing import List, Union, Tuple

import dask.bag as DB
from dask.distributed import progress

from fragmenstein import Victor
from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.steps.adapt_input import InputAdapter
from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.utils.io_utils import apply_func_to_files
from fragmenstein.utils.parallel_utils import get_parallel_client


class ErrorInComputation():
    pass


class NotComputedYet():
    pass


class CombineMerge_Base(ABC, InputAdapter):
    RESULT_PICKLE_TEMPLATE  = "%s.final.pickle"
    def __init__(self, output_path, template=None, templates_dir=None, template_pattern=None,
                 use_dask=False, verbose=False, *args, **kwargs):

        self.output_path = output_path
        self.template = template

        if self.template:
            assert templates_dir is None, "Error, if one template provided, templates_dir should be None"
            assert template_pattern is None, "Error, if one template provided, template_pattern should be None"
        else:
            assert templates_dir is not None, "Error, if no template provided, templates_dir should be provided instead"
            assert template_pattern is not None, "Error, if no template provided, template_pattern should be provided instead"

        self.template_pattern = template_pattern
        self.templates_dir = templates_dir
        self.use_dask = use_dask
        self.verbose = verbose

    @staticmethod
    def get_examples_init_params():
        return dict(
            output_path="./output_test_combineMerge",
            template=os.path.abspath(os.path.join(Xchem_info.examples_dir, "template.pdb"))
        )

    @staticmethod
    def get_examples_combine_params():
        data_dir = os.path.abspath(os.path.join(Xchem_info.examples_dir, "hit_mols"))
        fnames = [os.path.join(data_dir, "Mpro-x0678.mol"), os.path.join(data_dir, "Mpro-x0434.mol")]
        return dict(
            list_of_fragments=[[Compound.MolFromFile(fname, fname.split("-")[-1].split(".")[0]) for fname in fnames]]
        )

    @staticmethod
    def get_examples_place_params():
        data_dir = os.path.abspath(os.path.join(Xchem_info.examples_dir, "hit_mols"))
        fnames = [os.path.join(data_dir, "Mpro-x0107.mol"), os.path.join(data_dir, "Mpro-x0434.mol")]
        smi = "CC(NC(=O)CCl)c1cccc(Cl)c1"
        return dict(
            list_of_smi_fragments=[
                (smi, [Compound.MolFromFile(fname, fname.split("-")[-1].split(".")[0]) for fname in fnames])]
        )



    def get_final_results_name(self, merge_id, outdir=None):
        if not outdir:
            outdir = self.output_path
        return os.path.join(outdir, merge_id, CombineMerge_Base.RESULT_PICKLE_TEMPLATE%merge_id )

    def load_final_results(self, merge_id, outdir=None) -> Union[NotComputedYet, ErrorInComputation, Compound]:
        if not outdir:
            outdir = self.output_path
        try:
            with open(self.get_final_results_name(merge_id, outdir), "rb") as f:
                return pickle.load(f)
        except IOError:
            return NotComputedYet()

    def tryOneCombine(self, fragments: Union[List[str]]):  # add template as an optional parameter
        return self._tryOneGeneric(fragments)

    def tryOnePlace(self, smi: str, fragments: Union[List[str]]):  # add template as an optional parameter
        return self._tryOneGeneric(fragments, smi=smi)


    def getMergeId(self, fragment_ids, smi=None):
        merge_id = "-".join(fragment_ids)  # The order of the frag_ids is important
        if smi:
            merge_id = merge_id + "_" + smi  # It is important that smi goes after fragments for re.match in scoring

        merge_id = Victor.slugify(merge_id) # Very important since victor behaviour with names is differnt for combine and merge
        return merge_id

    def guessTemplate(self, fragment_ids):
        if self.template is None:  # If no global template, take the first fragment pdb as template

            template_fnames = apply_func_to_files(self.templates_dir, self.template_pattern, lambda x: x,
                                                  ids_to_check=fragment_ids)
            assert len(template_fnames)>0, "Error, template not found for fragment ids: %s"%str(fragment_ids)
            template_fnames = sorted(template_fnames)
            template = template_fnames[0]
            alternative_templates = template_fnames[1:]
            # template_fname = Victor.closest_hit(pdb_filenames=[f'{mpro_folder}/Mpro-{i}_0/Mpro-{i}_0_bound.pdb' for i in hit_codes], target_resi=145, target_chain='A',  target_atomname='SG', ligand_resn='LIG')
        else:
            template = self.template
            alternative_templates = []
        return template, alternative_templates

    def _tryOneGeneric(self, fragments, smi=None, *args, **kwargs) -> List[Compound] :

        fragments_dict = self.adapt_dict_or_compoundsList(fragments)

        frags, fragment_ids, bits_ids = [], [], []
        for frag in fragments_dict.values():
            frags.append(frag)
            fragment_ids.append(frag.primitiveId)
            bits_ids.append( frag.molId )

        merge_id = self.getMergeId( bits_ids, smi=smi)
        templateFname, alternative_fnames = self.guessTemplate(fragment_ids)

        prev_results = self.load_final_results(merge_id)
        if not isinstance(prev_results, NotComputedYet):
            return prev_results

        with tempfile.TemporaryDirectory() as tmp:
            wdir = os.path.join(tmp, merge_id)
            if not os.path.exists(wdir):
                os.mkdir(wdir)

            result_list = self.tryOneGeneric(merge_id, templateFname, frags, wdir, smi=smi)
            if not isinstance(result_list, ErrorInComputation) and len(result_list) == 0:
                result_list = ErrorInComputation()
            # print( result_list)
            # print(self.get_final_results_name(merge_id, outdir=Victor.work_path))
            # input(wdir + "  -> enter")

            with open(self.get_final_results_name(merge_id, outdir=wdir), "wb") as f:
                pickle.dump(result_list, f)

            dest_dir = os.path.join(self.output_path, merge_id)

            if os.path.isdir( dest_dir ):
                shutil.rmtree(dest_dir)
            shutil.copytree(os.path.join(wdir, merge_id), dest_dir )

        return result_list


    @abstractmethod
    def tryOneGeneric(self, merge_id, templateFname, fragments: List[Compound], wdir, smi: str = None, *args, **kwargs) -> List[Compound]:
        raise  NotImplementedError

    def applyGeneric(self, list_of_arguments, mapFunction):

        def keep_fun(x):
            return not isinstance(x, ErrorInComputation)

        if self.use_dask:
            dask_client = get_parallel_client()
            results_future = DB.from_sequence(list_of_arguments).map(mapFunction).filter(keep_fun)
            results = dask_client.compute(results_future)  # , scheduler='single-threaded')
            if self.verbose:
                progress(results)
            results = results.result()
        else:
            results = filter(keep_fun,
                                  map(mapFunction, list_of_arguments))
        results = chain.from_iterable( results)
        results = list( results )
        return results

    def applyCombine(self, list_of_fragments: List[List[Compound]]):
        '''

        :param list_of_fragments:
        :return:
        '''
        return self.applyGeneric(list_of_fragments, self.tryOneCombine)

    def applyPlace(self, list_of_smi_fragments: List[Tuple[str, List[Compound]]]):

        return self.applyGeneric(list_of_smi_fragments, lambda smi_frags: self.tryOnePlace(*smi_frags))
