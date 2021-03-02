import os
import sys
from rdkit import Chem


from collections import defaultdict

from fragmenstein.external.uploadToFragalysis.FragalysisFormater import FragalysisFormater
from fragmenstein.protocols.combineMerge_fragmensteinDefault import CombineMerge_FragmensteinDefault
from fragmenstein.protocols.hitsPreprocess_fragmentationBrics import HitsPreprocess_fragmentationBRICS
from fragmenstein.protocols.hitsPreprocess_permutations import HitsPreprocess_permutations
from fragmenstein.protocols.loadInput_XchemDefault import LoadInput_XchemDefault
from fragmenstein.protocols.score_combinedDefault import Score_CombinedDefault
from fragmenstein.utils.config_manager import ConfigManager

# output_dir = "~/oxford/tools/Fragmenstein/output"
# MAX_ATTEMPS = 2000
#
# template = "~/oxford/myProjects/diamondCovid/data/nsp13/aligned/nsp13-x0020_0B/nsp13-x0020_0B_apo-desolv.pdb"
# template_xchemId = "x0020"
#
# data_root_dir = "~/oxford/myProjects/diamondCovid/data/nsp13/aligned"
#
# #hit_ids = "x0041_0A,x0176_0A,x0176_1B,x0276_0A,x0309_0B,x0494_0A,x0494_0B".split(",") #"Site 7 Site 2 inner" WORKS
# #hit_ids = "x0058_0B,x0116_0B,x0306_0B,x0499_0B".split(",") #"Site 7 Site 2 outer" NOT WORKING
# #hit_ids = "x0020_0B,x0029_0B,x0257_0B".split(",") #"Site 2 B1 RNA-3'" NOT WORKING Flys away
# #hit_ids = "x0494_0A,x0494_0B,x0020_0B,x0029_0B,x0257_0B,x0041_0A,x0176_0A,x0176_1B,x0276_0A,x0309_0B,x0494_0A,x0494_0B".split(",") #"Site 7 all -> WORKS
# hit_ids = "x0176_0B,x0246_0B,x0438_0B".split(",") #""Site 1 - A Nucleoside" WORKS . Working on it
#
# sdf_outname = os.path.join(output_dir + "_site1_results.sdf")

RANDOM_SEED = 121
def main(data_root_dir, hit_ids, output_dir, preprocess_mode=None, template=None, template_xchemId=None, max_attemps=None, *args, **kwargs):

    data_root_dir = os.path.expanduser(data_root_dir)
    output_dir = os.path.expanduser(output_dir)

    if not os.path.exists(data_root_dir):
        os.mkdir(data_root_dir)
    sdf_outname = os.path.join(output_dir, ",".join(hit_ids)+".sdf")
    output_dir = os.path.join(output_dir, "fragmenstein" )

    loader = LoadInput_XchemDefault(data_root_dir, **LoadInput_XchemDefault.default_params_xchem() ) #TODO. Template should be found here
    fragments= list(filter( lambda frag: frag.molId in hit_ids, loader.fragments))

    proprocesser = None
    if preprocess_mode:
        print("preprocessing fragments")
        if preprocess_mode=="BRICS_decomposition":
            proprocesser = HitsPreprocess_fragmentationBRICS(fragments, random_seed= RANDOM_SEED)
            fragsCombin_iter =  proprocesser.yield_combinations(take_n_random=max_attemps)
        elif preprocess_mode=="SMARTS_decomposition":
            raise NotImplementedError()
        elif preprocess_mode == "permutations":
            proprocesser = HitsPreprocess_permutations(fragments, random_seed= RANDOM_SEED)
            fragsCombin_iter =  proprocesser.yield_combinations(max_num_elems=3, take_n_random=max_attemps)
        else:
            raise ValueError("Error, not implemented option preprocess_mode=%s"%preprocess_mode)
    else:
        fragsCombin_iter = [fragments]

    if template is None:
        id_template_iter = loader.find_templates( filter_ids=[hit_ids[0]])
        template_id, template = list( id_template_iter)[0]
        if not template_xchemId:
            template_xchemId = template_id
    else:
        template = os.path.expanduser(template)
        assert template_xchemId is not None, "Error, template_xchemId should not be None if template is not None "

    results = CombineMerge_FragmensteinDefault(output_path=output_dir, template=template, use_dask =ConfigManager.N_CPUS > 1).applyCombine(fragsCombin_iter)
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
        smi = Chem.MolToSmiles(generated_molecule, )
        #TODO: use MolHash
        if smi in already_available: continue
        already_available.add(smi)
        frags = score["fragments"]
        if proprocesser:
            frags = [ proprocesser.getOrinalFragmentId(frag) for frag in frags]
        generated_molecule.parents = frags
        proposed_mols[merge_id] = [generated_molecule, hit_ids]
        scores_dict[merge_id] = score

    scorer = Score_CombinedDefault(fragments_dir=data_root_dir, to_score_dir=output_dir,  selected_fragment_ids=hit_ids,
                                   ** Score_CombinedDefault.default_params_xchem())

    scores_dict = scorer.compute_scores(proposed_mols, already_computed_scores=scores_dict)

    # input( scores_dict )

    mols_list, metadata_list = zip(* scores_dict.values() )

    i = 0
    def get_simplified_mol_name(mol_id): #TODO: Move it within fragalysis??
        '''
        Reduce the mol_name for merges over fragments of fragments
        :param mol_id:
        :return:
        '''
        nonlocal i
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

        mol_id = "%d_"%i
        i+=1
        for fragId_chainId in sorted(ids_dict):
            mol_id += fragId_chainId
            for bit_id in sorted(ids_dict[fragId_chainId]):
                mol_id+="_"+",".join(bit_id)
        return mol_id

    for mol_name, (mol, score) in scores_dict.items():

        mol.SetProp("_Name", get_simplified_mol_name(mol_name))
        if proprocesser:
            score["fragments"] = sorted(set([ proprocesser.getOrinalFragmentId(frag) for frag in score["fragments"] ]))

    frag_writer = FragalysisFormater(ref_pdb_xchemId=template_xchemId)
    frag_writer.write_molsList_to_sdf(sdf_outname, mols_list, metadata_list)
    print("Results written to %s"%sdf_outname)

if __name__ == "__main__":
    from fragmenstein.utils.cmd_parser import ArgumentParser
    parser = ArgumentParser(prog="protocol_bricsFragmenstein", description="computes fragmenstein over BRICS fragmentation")
    parser.add_argument("-d", "--data_root_dir", type=str, help="The Xchem root dir for data, typically target_name/aligned/ ", required=True)
    parser.add_argument("-f", "--hit_ids", type=str, nargs="+", help="The hits ids to use  in the form x0020", required=True)
    parser.add_argument("-o", "--output_dir", type=str, help="The directory where results will be saved", required=True)
    parser.add_argument("-t", "--template", type=str, help="The path to a template pdb. If not provided, the first hit would be used", required=False)
    parser.add_argument("-x", "--template_xchemId", type=str, help="The xchem id that would be used for reference pdb in fragalysis", required=False)
    parser.add_argument("-p", "--preprocess_mode", type=str, choices=["BRICS_decomposition", "permutations", "SMARTS_decomposition"], default= None,
                        help="How to preprocess the fragments. XXX_decomposition breaks the fragments into bits.", required=False)
    parser.add_argument("-m", "--max_attemps", type=int, help="The number of maximun random attempts", required=False, default=None)

    args =vars( parser.parse_args())
    main( ** args)
    print("\nmain DONE!\n")

    '''

N_CPUS=1 python -m examples.protocols_bricsFragmenstein -d ~/oxford/myProjects/diamondCovid/data/nsp13/aligned -f x0176_0B x0246_0B x0438_0B -o ~/oxford/tools/Fragmenstein/output -m 10


python -m fragmenstein.external.condor_queue.send_to_condor --env_vars EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json DASK_WORKER_MEMORY=4GB --ncpus 8 "/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python -m examples.protocols_combineFragmenstein --n_cpus 8 -d /data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/aligned -o /data/xchem-fragalysis/sanchezg/oxford/tools/Fragmenstein/output_nsp13_Site7_mixed -f x0116_0B x0309_0B x4094_0B"

    '''