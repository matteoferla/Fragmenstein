import os
import sys
from rdkit import Chem


from collections import defaultdict

from fragmenstein.external.uploadToFragalysis.FragalysisFormater import FragalysisFormater
from fragmenstein.protocols.combineMerge_fragmensteinDefault import CombineMerge_FragmensteinDefault
from fragmenstein.protocols.fragmentation_brics import Fragmentator_BRICS
from fragmenstein.protocols.loadInput_XchemDefault import LoadInput_XchemDefault
from fragmenstein.protocols.score_combinedDefault import Score_CombinedDefault
from fragmenstein.utils.config_manager import ConfigManager

OUTPUT_PATH = "/home/ruben/oxford/tools/Fragmenstein/output"
MAX_ATTEMPS = 2000

template = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned/nsp13-x0020_0B/nsp13-x0020_0B_apo-desolv.pdb"
template_xchemId = "x0020"

hits_root_dir = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned"

hit_ids = [ "x0475_0A","x0509_0A" ] #"Site 20 & 22" WORKS but ugly
#hit_ids = "x0169_0B,x0290_0B,x0707_0B".split(",") #"Site 4-C1 RNA-5'" NOT WORKING
#hit_ids = "x0041_0A,x0176_0A,x0176_1B,x0276_0A,x0309_0B,x0494_0A,x0494_0B".split(",") #"Site 7 Site 2 inner" WORKS
#hit_ids = "x0058_0B,x0116_0B,x0306_0B,x0499_0B".split(",") #"Site 7 Site 2 outer" NOT WORKING
#hit_ids = "x0020_0B,x0029_0B,x0257_0B".split(",") #"Site 2 B1 RNA-3'" NOT WORKING
#hit_ids = "x0494_0A,x0494_0B,x0020_0B,x0029_0B,x0257_0B,x0041_0A,x0176_0A,x0176_1B,x0276_0A,x0309_0B,x0494_0A,x0494_0B".split(",") #"Site 7 all -> WORKS
# hit_ids = "x0176_0B,x0246_0B,x0438_0B".split(",") #""Site 1 - A Nucleoside" WORKS
# hit_ids = "x0020_0B,x0029_0B,x0257_0B".split(",") #""Site 1 - RNA-3' Site 1

sdf_outname = os.path.join(OUTPUT_PATH+"_site20&22_results.sdf")


def main(hit_ids = hit_ids):

    loader = LoadInput_XchemDefault( hits_root_dir, fragIds_pattern = LoadInput_XchemDefault.get_examples_init_params()["fragIds_pattern"]) #TODO. Template should be found here
    fragments= list(filter( lambda frag: frag.molId in hit_ids, loader.fragments))
    fragmentator = Fragmentator_BRICS( fragments , random_seed= 121)
    fragsCombin_iter =  fragmentator.yield_bits_combinations(take_n_random=MAX_ATTEMPS) #Refactored till here. Continue below

    results = CombineMerge_FragmensteinDefault(output_path=OUTPUT_PATH, template=template, use_dask = ConfigManager.N_CPUS>1).applyCombine( fragsCombin_iter)
    print("RESULTS:")
    print(results)

    if len(results) == 0:
      print("EXECUTION WAS NOT SUCCESSFUL.")
      sys.exit(0)

    proposed_mols={}
    scores_dict = {}
    already_available= set([])
    for merge_id, (score, generated_molecule) in results:
        #Heuristic filter for uniqueness
        smi = Chem.MolToSmiles(generated_molecule)
        if smi in already_available: continue
        already_available.add(smi)
        frags = score["fragments"]
        frags = [ fragmentator.getOrinalFragmentId(frag) for frag in frags]
        generated_molecule.parents = frags
        proposed_mols[merge_id] = [generated_molecule, hit_ids]
        scores_dict[merge_id] = score

    scorer = Score_CombinedDefault(fragments_dir=hits_root_dir, to_score_dir=OUTPUT_PATH, ** Score_CombinedDefault.default_params_xchem())
    scores_dict = scorer.compute_scores(proposed_mols, already_computed_scores=scores_dict)

    mols_list, metadata_list = zip(* scores_dict.values() )


    def get_simplified_mol_name(mol_id): #TODO: Move it within fragalysis??
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

    for mol_name, (mol, score) in scores_dict.items():

        mol.SetProp("_Name", get_simplified_mol_name(mol_name))
        score["fragments"] = [ fragmentator.getOrinalFragmentId(frag) for frag in score["fragments"] ]

    frag_writer = FragalysisFormater(ref_pdb_xchemId=template_xchemId)
    frag_writer.write_molsList_to_sdf("results_prueba.sdf", mols_list, metadata_list)




if __name__ == "__main__":
    main()
    print("\nmain DONE!\n")
    '''

N_CPUS=1 python -m examples.protocols_bricsFragmenstein

    '''