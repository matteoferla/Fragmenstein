import json
import logging
import os
import copy
from typing import List

from collections import defaultdict
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
                 merging_mode="none_permissive", verbose= None, save_pymol=False):

        if verbose is None:
            verbose = ConfigManager.VICTOR_VERBOSE

        self.save_pymol = save_pymol
        self.merging_mode = merging_mode
        super().__init__( output_path=output_path, template=template, templates_dir=templates_dir,
                          template_pattern=template_pattern, use_dask=use_dask, verbose=verbose)

        self._setVerbosity()

    def _setVerbosity(self):
        if self.verbose:
            Victor.enable_stdout(level=logging.DEBUG)
            Victor.quick_renanimation = True
        else:
            Victor.enable_stdout(level=logging.CRITICAL)
            Victor.enable_logfile("/dev/null", level=logging.CRITICAL)
            Victor.journal.setLevel(logging.CRITICAL)
            from rdkit_to_params import Params
            Params.log.setLevel("ERROR")
        Victor.error_to_catch = NotImplementedError

    def _resolveCovalentInfo(self, fragments, templateFname):
        covalent_info = None
        for frag in fragments:
            if frag.covalent_info is not None:
                covalent_info = frag.covalent_info #Pick one covalently bound residue to define restrains-

        if covalent_info is None:
            distManager = PdbDistanceManager(templateFname, limit_to_resname="CYS")
            for frag in fragments:
                chainId_resId_resname = distManager.find_closest_residue(frag)
                if chainId_resId_resname:
                    covalent_info = {'covalent_resi': "".join(reversed(chainId_resId_resname[:2])), 'covalent_resn':'CYS' }
                    break
        return covalent_info

    def tryOneGeneric(self, merge_id, templateFname, fragments: List[Compound], wdir, smi: str = None,
                      alternative_templateFnames=[], *args, **kwargs):

        ExternalToolImporter.import_tool("pyrosetta", ["pyrosetta"])

        covalent_info = self._resolveCovalentInfo( fragments, templateFname)

        Victor.work_path = wdir

        minimized_mol = None
        unminimized_mol_pdbblock = None
        v = Victor(hits=fragments, pdb_filename=templateFname, **covalent_info)
        try:
            if smi is None:
                v.combine( long_name=merge_id ) #Warning long_name= merge_id will mutate merge_id to replace "_" -> "-", so should be obtained from self.getMergeId(fragIds)
            else:
                assert  Chem.MolFromSmiles(smi) is not None, "Error, invalid smiles: %s"%smi
                v.place( smi, long_name=merge_id, merging_mode=self.merging_mode )
            minimized_mol = copy.deepcopy(v.minimised_mol)
            _unminimized_mol = copy.deepcopy(v.monster.positioned_mol)
            Chem.SanitizeMol(minimized_mol)
            Chem.SanitizeMol(_unminimized_mol)
            unminimized_mol_pdbblock = Chem.MolToPDBBlock(_unminimized_mol)

        except Exception as e:
            print(e)

        if minimized_mol is not None and unminimized_mol_pdbblock is not None:
            if self.save_pymol:
                v.make_pse()

            metadata_dict = v.summarise()
            metadata_dict = _FragmensteinScorer.old_scoring_fun(metadata_dict)[-1]
            metadata_dict = _FragmensteinScorer.new_scoring_fun(metadata_dict)[-1]

            ori_frag_ids = set([ x.primitiveId for x in fragments])

            derived_frags = metadata_dict["regarded"]
            regarded_fragIds = sorted(set([ oriFragId for oriFragId in ori_frag_ids
                                                    if any([oriFragId in devFrag
                                                            for devFrag in derived_frags])]))

            atomic_mapping = self.find_atoms_mapping(_unminimized_mol, filter(lambda x: x.primitiveId in regarded_fragIds, fragments))

            metadata_dict["atomic_mapping"] = atomic_mapping
            metadata_dict["fragments"] = regarded_fragIds
            metadata_dict["ref_pdb"] = templateFname
            minimized_mol = Compound( minimized_mol, molId=merge_id, parents= fragments)
            minimized_mol.ref_pdb = templateFname
            minimized_mol.metadata = metadata_dict
            minimized_mol.ref_molIds =  regarded_fragIds
            minimized_mol.atomic_mapping = atomic_mapping
            result = [ minimized_mol ]

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