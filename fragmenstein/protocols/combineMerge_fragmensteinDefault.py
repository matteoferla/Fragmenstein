import json
import os
import pickle
import tempfile

import shutil
from typing import List, Dict, Union, Tuple

import logging
from rdkit import Chem
import dask.bag as DB
from dask.distributed import progress

from fragmenstein import Victor
from fragmenstein.external import ExternalToolImporter
from fragmenstein.protocols.adapt_input import InputAdapter
from fragmenstein.scoring._fragmenstein_scoring import _FragmensteinScorer
from fragmenstein.utils.compound import Compound
from fragmenstein.utils.config_manager import ConfigManager
from fragmenstein.utils.io_utils import apply_func_to_files
from fragmenstein.utils.parallel_utils import get_parallel_client
from fragmenstein.utils.pdb_utils import PdbDistanceManager


class ErrorInComputation():
    pass

class NotComputedYet():
    pass

class CombineMerge_FragmensteinDefault(InputAdapter):

    RESULT_METADATA_PATTERN= "%s.scores.json"
    def __init__(self, output_path, template=None, templates_dir=None, template_pattern=None, merging_mode="none_permissive", use_dask=False,
                 verbose= ConfigManager.VICTOR_VERBOSE, save_pymol=False):

        self.save_pymol = save_pymol
        self.verbose = verbose
        self.output_path = output_path
        self.template = template
        self.merging_mode = merging_mode

        if self.template:
            assert templates_dir is None, "Error, if one template provided, templates_dir should be None"
            assert template_pattern is None, "Error, if one template provided, template_pattern should be None"
        else:
            assert templates_dir is not None, "Error, if no template provided, templates_dir should be provided instead"
            assert template_pattern is not None, "Error, if no template provided, template_pattern should be provided instead"

        self.template_pattern = template_pattern
        self.templates_dir = templates_dir
        self.use_dask = use_dask


    @staticmethod
    def get_examples_init_params():
        return dict(
                    output_path="./output_test_combineMerge",
                    template= os.path.abspath( os.path.join(__file__, "../../mpro/data/template.pdb") )
        )

    @staticmethod
    def get_examples_combine_params():
        data_dir = os.path.abspath( os.path.join(__file__, "../../mpro/data/hit_mols"))
        fnames = [ os.path.join(data_dir, "Mpro-x0678.mol"), os.path.join(data_dir, "Mpro-x0434.mol")  ]
        return dict(
            list_of_fragments= [ [Compound.MolFromFile(fname, fname.split("-")[-1].split(".")[0]) for fname in fnames ] ]
        )

    @staticmethod
    def get_examples_place_params():
        data_dir = os.path.abspath( os.path.join(__file__, "../../mpro/data/hit_mols"))
        fnames = [ os.path.join(data_dir, "Mpro-x0107.mol"), os.path.join(data_dir, "Mpro-x0434.mol")  ]
        smi = "CC(NC(=O)CCl)c1cccc(Cl)c1"
        return dict(
            list_of_smi_fragments= [ (smi, [Compound.MolFromFile(fname, fname.split("-")[-1].split(".")[0]) for fname in fnames ]) ]
        )

    def get_final_results_name(self, merge_id, outdir=None):
        if not outdir:
            outdir = self.output_path
        return os.path.join(outdir, merge_id, merge_id + ".final.pickle")

    def load_final_results(self, merge_id, outdir=None) -> Union[NotComputedYet, ErrorInComputation, Compound]:
        if not outdir:
            outdir = self.output_path
        try:
            with open(self.get_final_results_name(merge_id, outdir), "rb") as f:
                return pickle.load(f)
        except IOError:
            return NotComputedYet()

    def tryOneCombine(self, fragments:  Union[List[str]]): #add template as an optional parameter
        return self.tryOneGeneric(fragments)

    def tryOnePlace(self, smi: str, fragments:  Union[List[str]]): #add template as an optional parameter
        return self.tryOneGeneric(fragments, smi= smi)

    def tryOneGeneric(self, fragments:  Union[List[str]], smi: str= None): #add template as an optional parameter

        ExternalToolImporter.import_tool("pyrosetta", ["pyrosetta"])
        fragments_dict = self.adapt_dict_or_compoundsList( fragments )

        frags, fragment_ids = [], []
        covalent_info = None
        for fragId, frag in fragments_dict.items():
            frag.SetProp('_Name', fragId)
            frags.append(frag)
            fragment_ids.append( fragId )
            if frag.covalent_info is not None:
                covalent_info = frag.covalent_info #Pick one covalently bound residue to define restrains-

        merge_id = "-".join(fragment_ids) #The order of the frag_ids is important
        if smi:
            merge_id =  merge_id +"_"+ smi #It is important that smi goes after fragments for re.match in scoring

        merge_id = Victor.slugify(merge_id) #Very important since victor behaviour with names is differnt for combine and merge


        if self.template is None: #If no global template, take the first fragment pdb as template

            template_fnames = apply_func_to_files(self.templates_dir, self.template_pattern, lambda x: x, ids_to_check=fragment_ids)
            # print(self.templates_dir, self.template_pattern, fragment_ids)
            template = sorted( template_fnames)[0]

            # template_fname = Victor.closest_hit(pdb_filenames=[f'{mpro_folder}/Mpro-{cur_mol_num}_0/Mpro-{cur_mol_num}_0_bound.pdb' for cur_mol_num in hit_codes], target_resi=145, target_chain='A',  target_atomname='SG', ligand_resn='LIG')
        else:
            template = self.template


        if covalent_info is None:
            distManager = PdbDistanceManager(template, limit_to_resname="CYS")
            for frag in frags:
                chainId_resId = distManager.find_closest_residue(frag)
                if chainId_resId:
                    covalent_info = {'covalent_resi': "".join(reversed(chainId_resId)), 'covalent_resn':'CYS' }
                    break

        prev_results = self.load_final_results(merge_id)
        if not isinstance(prev_results, NotComputedYet ):
            return prev_results
        #else: compute
    
        with tempfile.TemporaryDirectory() as tmp:
            Victor.work_path = tmp
            Victor.error_to_catch = NotImplementedError
            if self.verbose:
                Victor.enable_stdout(level=logging.DEBUG)
                Victor.enable_logfile(level=logging.DEBUG)
                # Victor.error_to_catch = NotImplementedError
            else:
                Victor.enable_stdout(level=logging.CRITICAL)
                Victor.enable_logfile("/dev/null", level=logging.CRITICAL)
                Victor.journal.setLevel(logging.CRITICAL)
                # Victor.error_to_catch = NotImplementedError

            generated_molecule = None
            v = Victor(hits=frags, pdb_filename=template, **covalent_info)
            try:
                if smi is None:
                    v.combine( long_name=merge_id ) #Warning long_name= merge_id will mutate merge_id to replace "_" -> "-", so already done in prev line
                else:
                    assert  Chem.MolFromSmiles(smi) is not None, "Error, invalid smiles: %s"%smi
                    v.place( smi, long_name=merge_id, merging_mode=self.merging_mode )
                generated_molecule = v.minimised_mol
            except Exception as e:
                print(e)
                pass
    
            if generated_molecule is not None:
                if self.save_pymol:
                    v.make_pse()

                metadata_dict = v.summarise()
                metadata_dict = _FragmensteinScorer.old_scoring_fun(metadata_dict)[-1]
                metadata_dict = _FragmensteinScorer.new_scoring_fun(metadata_dict)[-1]

                metadata_dict["fragments"] = metadata_dict["regarded"]
                metadata_dict["ref_pdb"] = template

                generated_molecule = Compound( generated_molecule, molId =merge_id, parents= frags)
                generated_molecule.ref_pdb = template
                generated_molecule.metadata = metadata_dict
                generated_molecule.ref_molIds =  metadata_dict["regarded"]

                result = generated_molecule

            else:
                metadata_dict = {"error": v.error_msg }
                result = ErrorInComputation()

            if not os.path.exists(os.path.join(Victor.work_path, merge_id)):
                os.mkdir(os.path.join(Victor.work_path, merge_id))

            scores_json_basename = CombineMerge_FragmensteinDefault.RESULT_METADATA_PATTERN % merge_id
            with open(os.path.join(Victor.work_path, merge_id, scores_json_basename), "w") as f:
                json.dump(metadata_dict, f)
            with open( self.get_final_results_name(merge_id, outdir=tmp), "wb") as f:
                pickle.dump(result, f)

            shutil.copytree(os.path.join(Victor.work_path, merge_id), os.path.join(self.output_path, merge_id))
        return result


    def applyGeneric(self, list_of_arguments, mapFunction ):

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
            results = list(filter(keep_fun,
                                  map( mapFunction, list_of_arguments)))
        return results

    def applyCombine(self, list_of_fragments: List[ List[Compound] ]):
        '''

        :param list_of_fragments:
        :return:
        '''
        return self.applyGeneric(list_of_fragments, self.tryOneCombine)

    def applyPlace(self, list_of_smi_fragments: List[ Tuple[str, List[Compound]]] ):

        return self.applyGeneric(list_of_smi_fragments, lambda smi_frags: self.tryOnePlace(*smi_frags))

def test_applyCombine():
    placer = CombineMerge_FragmensteinDefault( **CombineMerge_FragmensteinDefault.get_examples_init_params(), use_dask=True)
    results = placer.applyCombine( **CombineMerge_FragmensteinDefault.get_examples_combine_params())
    print("RESULTS applyCombine:")
    print( results)


def test_applyPlace():
    placer = CombineMerge_FragmensteinDefault( **CombineMerge_FragmensteinDefault.get_examples_init_params(), use_dask=False)
    results = placer.applyPlace( **CombineMerge_FragmensteinDefault.get_examples_place_params())
    print("RESULTS applyPlace:")
    print( results)

if __name__ == "__main__":

    test_applyCombine()
    test_applyPlace()

    '''
python -m fragmenstein/protocols/combineMerge_fragmensteinDefault.py
    '''