import json
import logging
import os
from typing import List

from rdkit import Chem

from fragmenstein import Victor
from fragmenstein.external import ExternalToolImporter
from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.steps.combineMerge_abstract import ErrorInComputation, CombineMerge_Base
from fragmenstein.scoring._fragmenstein_scoring import _FragmensteinScorer
from fragmenstein.utils.config_manager import ConfigManager
from fragmenstein.utils.pdb_utils import PdbDistanceManager


class CombineMerge_FragmensteinDefault( CombineMerge_Base ):

    RESULT_METADATA_PATTERN= "%s.scores.json"

    def __init__(self, output_path, template=None, templates_dir=None, template_pattern=None, use_dask=False,
                 merging_mode="none_permissive", verbose= ConfigManager.VICTOR_VERBOSE, save_pymol=False):

        self.save_pymol = save_pymol
        self.merging_mode = merging_mode
        super().__init__( output_path=output_path, template=template, templates_dir=templates_dir,
                          template_pattern=template_pattern, use_dask=use_dask, verbose=verbose)

    def tryOneGeneric(self, merge_id, templateFname, fragments: List[Compound], wdir, smi: str = None, *args, **kwargs):

        ExternalToolImporter.import_tool("pyrosetta", ["pyrosetta"])

        covalent_info = None
        for frag in fragments:
            if frag.covalent_info is not None:
                covalent_info = frag.covalent_info #Pick one covalently bound residue to define restrains-

        if covalent_info is None:
            distManager = PdbDistanceManager(templateFname, limit_to_resname="CYS")
            for frag in fragments:
                chainId_resId = distManager.find_closest_residue(frag)
                if chainId_resId:
                    covalent_info = {'covalent_resi': "".join(reversed(chainId_resId)), 'covalent_resn':'CYS' }
                    break
    
        Victor.work_path = wdir
        if self.verbose:
            Victor.enable_stdout(level=logging.DEBUG)
        else:
            Victor.enable_stdout(level=logging.CRITICAL)
            Victor.enable_logfile("/dev/null", level=logging.CRITICAL)
            Victor.journal.setLevel(logging.CRITICAL)

        # Victor.error_to_catch = NotImplementedError

        generated_molecule = None
        v = Victor(hits=fragments, pdb_filename=templateFname, **covalent_info)
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
            metadata_dict["ref_pdb"] = templateFname

            generated_molecule = Compound( generated_molecule, molId=merge_id, parents= fragments)
            generated_molecule.ref_pdb = templateFname
            generated_molecule.metadata = metadata_dict
            generated_molecule.ref_molIds =  metadata_dict["regarded"]

            result = [ generated_molecule ]

        else:
            metadata_dict = {"error": v.error_msg }
            result = ErrorInComputation()

        if not os.path.exists(os.path.join(Victor.work_path, merge_id)):
            os.mkdir(os.path.join(Victor.work_path, merge_id))

        scores_json_basename = CombineMerge_FragmensteinDefault.RESULT_METADATA_PATTERN % merge_id
        with open(os.path.join(Victor.work_path, merge_id, scores_json_basename), "w") as f:
            json.dump(metadata_dict, f)

        return result

def test_applyCombine():
    combiner = CombineMerge_FragmensteinDefault( **CombineMerge_FragmensteinDefault.get_examples_init_params(), use_dask=True)
    results = combiner.applyCombine( **CombineMerge_FragmensteinDefault.get_examples_combine_params())
    print("RESULTS applyCombine:")
    print( results)


def test_applyPlace():
    placer = CombineMerge_FragmensteinDefault( **CombineMerge_FragmensteinDefault.get_examples_init_params(), use_dask=False)
    results = placer.applyPlace( **CombineMerge_FragmensteinDefault.get_examples_place_params())
    print("RESULTS applyPlace:")
    print( results)

if __name__ == "__main__":

    print("trying combine")
    test_applyCombine()
    print("trying place")
    test_applyPlace()

    '''
python -m fragmenstein.protocols.steps.combineMerge_fragmensteinDefault
    '''