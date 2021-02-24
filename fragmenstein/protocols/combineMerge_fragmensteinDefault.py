import json
import os
import tempfile

import shutil
from typing import List, Dict, Union, Tuple

from rdkit import Chem
import dask.bag as DB
from dask.distributed import progress

from fragmenstein import Victor
from fragmenstein.external import ExternalToolImporter
from fragmenstein.protocols.adapt_input import InputAddapter
from fragmenstein.scoring._fragmenstein_scoring import _FragmensteinScorer
from fragmenstein.utils.compound import Compound
from fragmenstein.utils.parallel_utils import get_parallel_client
from fragmenstein.utils.pdb_utils import PdbDistanceManager


class CombineMerge_FragmensteinDefault(InputAddapter):

    RESULT_METADATA_PATTERN= "%s.scores.json"
    def __init__(self, output_path, template=None, use_dask=False, verbose=False):


        self.verbose = verbose
        self.output_path = output_path
        self.template = template #TODO: if template is null, find template within
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
            merge_id = smi + merge_id

        merge_id = Victor.slugify(merge_id) #Very important since victor behaviour with names is differnt for combine and merge


        # if self.template is None:        #TODO: if template == None select the best one. Better before.
        #     template_fname = Victor.closest_hit(pdb_filenames=[f'{mpro_folder}/Mpro-{i}_0/Mpro-{i}_0_bound.pdb' for i in hit_codes],
        #                         target_resi=145,
        #                         target_chain='A',
        #                         target_atomname='SG',
        #                         ligand_resn='LIG')



        if covalent_info is None:
            distManager = PdbDistanceManager(self.template, limit_to_resname="CYS")
            for frag in frags:
                chainId_resId = distManager.find_closest_residue(frag)
                if chainId_resId:
                    covalent_info = {'covalent_resi': "".join(reversed(chainId_resId)), 'covalent_resn':'CYS' }
                    break
    
        scores_json_basename = CombineMerge_FragmensteinDefault.RESULT_METADATA_PATTERN % merge_id
        scores_json_fname = os.path.join(self.output_path, merge_id, scores_json_basename)

        # print( "\n"+scores_json_fname+ "\n")

        if os.path.exists(scores_json_fname):
            with open(scores_json_fname) as f:
                metadata_dict = json.load(f)

            mol_fname = os.path.join(self.output_path, merge_id, merge_id + ".minimised.mol")
            if os.path.exists(mol_fname):
                generated_molecule = Chem.MolFromMolFile(mol_fname)
                generated_molecule = Compound( generated_molecule, molId=merge_id, parents= frags)
                return merge_id, (metadata_dict, generated_molecule)
            else:
                return None
    
        with tempfile.TemporaryDirectory() as tmp:
            Victor.work_path = tmp

            if self.verbose:
                import logging
                Victor.enable_stdout(level=logging.DEBUG)
                Victor.enable_logfile(level=logging.DEBUG)
                # Victor.error_to_catch = NotImplementedError

            generated_molecule = None
            v = Victor(hits=frags, pdb_filename=self.template,
                       **covalent_info)
            try:
                if smi is None:
                    v.combine( long_name=merge_id ) #Warning long_name= merge_id will mutate merge_id to replace "_" -> "-"
                else:
                    assert  Chem.MolFromSmiles(smi) is not None, "Error, invalid smiles: %s"%smi
                    v.place( smi, long_name=merge_id )
                generated_molecule = v.minimised_mol
            except Exception as e:
                print(e)
                pass
    
            if generated_molecule is not None:
                v.make_pse()
                metadata_dict = v.summarise()
                metadata_dict = _FragmensteinScorer.old_scoring_fun(metadata_dict)[-1]
                metadata_dict = _FragmensteinScorer.new_scoring_fun(metadata_dict)[-1]

                metadata_dict["fragments"] = fragment_ids
                generated_molecule = Compound( generated_molecule, molId =merge_id, parents= frags)
                result = merge_id, ( metadata_dict, generated_molecule)
            else:
                metadata_dict = None
                result = None
            if not os.path.exists(os.path.join(Victor.work_path, merge_id)):
                os.mkdir(os.path.join(Victor.work_path, merge_id))

            with open(os.path.join(Victor.work_path, merge_id, scores_json_basename), "w") as f:
                json.dump(metadata_dict, f)
            shutil.copytree(os.path.join(Victor.work_path, merge_id), os.path.join(self.output_path, merge_id))
    
        return result


    def applyGeneric(self, list_of_arguments, mapFunction ):

        if self.use_dask:
            dask_client = get_parallel_client()
            results_future = DB.from_sequence(list_of_arguments).map(mapFunction
                                                              ).filter(lambda x: x is not None)
            results = dask_client.compute(results_future)  # , scheduler='single-threaded')
            progress(results)
            results = results.result()
        else:
            results = list(filter(None.__ne__,
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
    placer = CombineMerge_FragmensteinDefault( **CombineMerge_FragmensteinDefault.get_examples_init_params(), initiaze_parallel_execution=True)
    results = placer.applyCombine( **CombineMerge_FragmensteinDefault.get_examples_combine_params())
    print("RESULTS applyCombine:")
    print( results)


def test_applyPlace():
    placer = CombineMerge_FragmensteinDefault( **CombineMerge_FragmensteinDefault.get_examples_init_params(), initiaze_parallel_execution=False)
    results = placer.applyPlace( **CombineMerge_FragmensteinDefault.get_examples_place_params())
    print("RESULTS applyPlace:")
    print( results)

if __name__ == "__main__":

    test_applyCombine()
    test_applyPlace()

    '''
python -m fragmenstein/protocols/combineMerge_fragmensteinDefault.py
    '''