import json
import random
import sys, os
import tempfile
import logging
import shutil

from itertools import chain
from itertools import combinations
from collections import defaultdict
from joblib import Parallel,  delayed
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

from fragmenstein import Victor

from fragmenstein.external.uploadToFragalysis.FragalysisFormater import FragalysisFormater
from fragmenstein.scoring._fragmenstein_scoring import _FragmensteinScorer
from fragmenstein.scoring.combined_scorer import CombineScorer
from fragmenstein.scoring.interactionBasedScorer import InteractionBasedScorer
from fragmenstein.scoring.propertiesScorer import PropertiesScorer
from fragmenstein.scoring.sucos import SuCOSComputer
from fragmenstein.scoring.xcos import XcosComputer
from fragmenstein.utils.config_manager import ConfigManager
from fragmenstein.utils.fragmentation_utils import split_mol_to_brics_bits

Victor.enable_stdout(level=logging.DEBUG)
# Victor.error_to_catch = NotImplementedError

import pyrosetta
pyrosetta.init(extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')


OUTPUT_PATH = "/home/ruben/oxford/tools/Fragmenstein/output"
MAX_ATTEMPS = 10
INTERACTIVE_RESULTS = False
N_JOBs = ConfigManager.N_CPUS

template = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned/nsp13-x0020_0B/nsp13-x0020_0B_apo-desolv.pdb"
template_xchemId = "x0020"

hits_root_dir = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned"

hits_pattern = hits_root_dir + "/nsp13-%(hit_id)s/nsp13-%(hit_id)s.mol"


# hit_ids = [ "x0475_0A","x0509_0A" ] #"Site 20 & 22" WORKS but ugly
#hit_ids = "x0169_0B,x0290_0B,x0707_0B".split(",") #"Site 4-C1 RNA-5'" NOT WORKING
#hit_ids = "x0041_0A,x0176_0A,x0176_1B,x0276_0A,x0309_0B,x0494_0A,x0494_0B".split(",") #"Site 7 Site 2 inner" WORKS
#hit_ids = "x0058_0B,x0116_0B,x0306_0B,x0499_0B".split(",") #"Site 7 Site 2 outer" NOT WORKING
#hit_ids = "x0020_0B,x0029_0B,x0257_0B".split(",") #"Site 2 B1 RNA-3'" NOT WORKING
#hit_ids = "x0494_0A,x0494_0B,x0020_0B,x0029_0B,x0257_0B,x0041_0A,x0176_0A,x0176_1B,x0276_0A,x0309_0B,x0494_0A,x0494_0B".split(",") #"Site 7 all -> WORKS
# hit_ids = "x0176_0B,x0246_0B,x0438_0B".split(",") #""Site 1 - A Nucleoside" WORKS
hit_ids = "x0020_0B,x0029_0B,x0257_0B".split(",") #""Site 1 - RNA-3' Site 1


def load_hit(hit_id):
  mol = Chem.MolFromMolFile( hits_pattern%{"hit_id":hit_id} )
  if mol.GetProp('_Name') == '':
      mol.SetProp('_Name', hit_id )
  return hit_id, mol


hitId_hit_list = list(map(load_hit, hit_ids))
bitId_bit_list = []
for hit_id, hit in hitId_hit_list:
    bits = split_mol_to_brics_bits(hit)
    bits = sorted( bits, key=lambda  x: ExactMolWt(x) )
    for i, bit in enumerate( bits):
        bit_id = hit_id + "_%d" % i
        bit = Chem.DeleteSubstructs(bit, Chem.MolFromSmiles('*'))  # TODO: if covalently interacting with the protein, dummy atoms before fragmentation should not be removed
        bitId_bit_list.append((bit_id, bit))


def powerset(iterable, min_num_frag=2, include_full=False):
    s = list(iterable)
    last = len(s)+1 if include_full else len(s)
    return list(chain.from_iterable(combinations(s, r) for r in range(min_num_frag, last)))

def draw_mols(mols):
    from matplotlib import pyplot as plt
    from rdkit.Chem import Draw
    fig, ax = plt.subplots(2,1)
    ax[0].imshow(Draw.MolsToGridImage(mols, molsPerRow=3))
    ax[1].imshow(Draw.MolsToGridImage([ Chem.MolFromSmiles(Chem.MolToSmiles(mol)) for mol in mols], molsPerRow=3))
    plt.show()

bitId_bit_list = powerset(bitId_bit_list)
random.shuffle(bitId_bit_list)
bitId_bit_list = [hitId_hit_list] + bitId_bit_list


def tryOneMerge(trial_hitId_hit_list):
    hit_ids, hits  = zip(*trial_hitId_hit_list)
    merge_id = "-" .join(hit_ids)

    for hit_id, hit in trial_hitId_hit_list: #Set prop should be done within subprocess ( delayed call) otherwise will be lost
        hit.SetProp('_Name', hit_id )

    scores_json_basename = merge_id+".scores.json"
    scores_json_fname = os.path.join(OUTPUT_PATH, merge_id, scores_json_basename)
    # scores_json_fname = os.path.join(OUTPUT_PATH, merge_id, "scores.json") #TODO: replace with the above expresion

    if os.path.exists(scores_json_fname):
        with open(scores_json_fname) as f:
            score = json.load(f)
        mol_fname = os.path.join(OUTPUT_PATH, merge_id, merge_id + ".minimized.mol")
        if  os.path.exists(mol_fname):
            generated_molecule = Chem.MolFromMolFile(mol_fname)
            return merge_id, (score, generated_molecule)
        else:
            return  None

    with tempfile.TemporaryDirectory() as tmp:
        Victor.work_path= tmp
        v = Victor.combine(hits=hits, pdb_filename=template,
                           # a random residue is still required for the constraint ref atom.
                           covalent_resi='8B', covalent_resn='CYS')
        generated_molecule = v.minimised_mol

        if generated_molecule is not None:
            # draw_mols( hitId_hit_list+ [ generated_molecule ] )
            v.make_pse()
            score = v.summarise()
            score = _FragmensteinScorer.old_scoring_fun(score)[-1]
            score = _FragmensteinScorer.new_scoring_fun(score)[-1]
            # fragments = sorted(set( [ "_".join(elem.split("_")[:2])[:-1] for elem in merge_id.split("-") ] ))
            fragments = sorted(set( [ "_".join(elem.split("_")[:2]) for elem in merge_id.split("-") ] ))

            score["fragments"] = fragments
            result = merge_id, (score, generated_molecule)  #TODO: old_scoring_fun and new_scoring_fun are not being saved
        else:
            score = None
            result = None
            if not os.path.exists(  os.path.join(Victor.work_path, merge_id) ):
                os.mkdir( os.path.join(Victor.work_path, merge_id) )
        with open( os.path.join(Victor.work_path, merge_id, scores_json_basename), "w" ) as f:
            json.dump( score, f)
        shutil.copytree(os.path.join(Victor.work_path, merge_id), os.path.join(OUTPUT_PATH, merge_id))

        return result

results = Parallel(n_jobs=N_JOBs, batch_size=1, backend="multiprocessing")(delayed(tryOneMerge)(bitId_bit_list[i]) for i in range(0, min(len(bitId_bit_list), MAX_ATTEMPS)))
results = list(filter(None.__ne__, results ) )


if len(results) == 0:
  print("EXECUTION WAS NOT SUCCESSFUL.")
  sys.exit(0)

proposed_mols={}
for merge_id, (score, generated_molecule) in results:
    proposed_mols[merge_id] = [generated_molecule, hit_ids]

def get_simplified_mol_name(mol_id):
    '''
    Reduce the mol_name for merges over fragments of fragments
    :param mol_id:
    :return:
    '''
    pieces = mol_id.split("-")
    ids_dict = defaultdict( set)
    for frag_id in pieces:
        split_fragId = frag_id.split("_")
        fragId_chainId = "-".join(split_fragId[:2])
        if len(split_fragId)==3:
            bit_id = split_fragId[2]
        else:
            bit_id = ""

        ids_dict[fragId_chainId].add(bit_id)

    mol_id = ""
    for fragId_chainId in sorted(ids_dict):
        mol_id += fragId_chainId
        for bit_id in sorted(ids_dict[fragId_chainId]):
            mol_id+="_"+",".join(bit_id)
    return mol_id

with tempfile.TemporaryDirectory() as tmp:
    fragment_id_pattern = r".*-(\w+)\.mol$"
    scorer1 = SuCOSComputer(fragments_dir=hits_root_dir, fragment_id_pattern=fragment_id_pattern, working_dir=tmp )
    scorer2 = PropertiesScorer(working_dir=tmp)
    scorer3 = XcosComputer(fragments_dir=hits_root_dir, fragment_id_pattern=fragment_id_pattern, working_dir=tmp )
    scorer4 = InteractionBasedScorer( fragments_dir=hits_root_dir, fragment_id_pattern="(.+)_bound\.pdb$",
                                      boundPdbs_to_score_dir=OUTPUT_PATH, boundPdbs_to_score_pattern= ".*-(\w+)_bound\.pdb$", working_dir=tmp )
    scorers_list = [scorer1, scorer2, scorer3, scorer4]
    scores_list = CombineScorer.computeScoreForMolecules(proposed_mols , scorers_objects_list=scorers_list, working_dir=tmp)

metadata_list = []
mols_list = []
for (mol_name, (score, mol)), additional_score in zip(results, scores_list):
    mol.SetProp("_Name", get_simplified_mol_name(mol_name))
    mol.SetProp("original SMILES", Chem.MolToSmiles(mol))
    mol.SetProp("ref_mols", ",".join(score["fragments"]))

    mols_list.append(mol)
    score.update(additional_score)
    metadata_list.append(score)


frag_writer = FragalysisFormater(ref_pdb_xchemId=template_xchemId)
frag_writer.write_molsList_to_sdf("results_prueba.sdf", mols_list, metadata_list)