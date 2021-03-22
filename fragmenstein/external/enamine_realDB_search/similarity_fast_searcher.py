import datetime
import json
import os
import sqlite3
import time
#os.environ["NUMBA_NUM_THREADS"]=4
import numba
from numba.types import Array
from collections import OrderedDict

import dask.bag as db
import numpy as np

from fragmenstein.external.enamine_realDB_search.common import get_fingerPrint_as_npBool, \
    decompressFingerprint_npStr
from fragmenstein.utils.cmd_parser import ArgumentParser
from fragmenstein.utils.parallel_utils import get_parallel_client



def jaccard_vectorized(x,y):
  intersections_count= x.astype(np.float32) @ y.T
  counts_x = np.sum(x, axis=-1)
  counts_y = np.sum(y, axis=-1)
  sums = counts_x.reshape(-1,1)+counts_y.reshape(-1,1).T
  union = sums - intersections_count
  return intersections_count / union


@numba.jit( nopython=True, cache=True)
def compute_2FingerPrints_similarity(query_fp, db_fp):

    n11 = np.sum((query_fp == 1) & (db_fp == 1))
    n00 = np.sum((query_fp == 0) & (db_fp == 0))
    jac = n11 / (query_fp.shape[0] - n00)
    return jac

numba_decompressFingerprint_npStr = numba.jit(decompressFingerprint_npStr, nopython=True, cache=True)

@numba.jit( nopython=True, cache=True)
def compute_FingerPrints_similarity_from_allbools(query_fp_matrix, db_many_fps_bool, file_num, n_hits_per_smi):

    n_query_fps = query_fp_matrix.shape[0]
    fp_bipts = query_fp_matrix.shape[1]
    matched_similarities = np.ones((n_query_fps, n_hits_per_smi)) * -1  # query_id, hit_num, similarity
    matched_ids = np.ones((n_query_fps, n_hits_per_smi, 2),
                          dtype=np.int64) * -1  # query_id, hit_num, [ file_id, hit_id]


    n_byes_in_db_fps = len(db_many_fps_bool)
    assert n_byes_in_db_fps >= fp_bipts, "Error, fingerprints_byte string does not contain fingerprints"
    assert n_byes_in_db_fps % fp_bipts ==0, "Error, somo corruption found in db_many_fps_bool"
    n_mols = n_byes_in_db_fps//fp_bipts

    n_mol_processed =0

    print_step = 10**int(np.ceil(np.log10((int(1e5) //n_query_fps))))
    for j in range(0, n_byes_in_db_fps, fp_bipts):
        bin_finPrint = db_many_fps_bool[j:j + fp_bipts]
        n_mol_processed += 1
        if n_mol_processed % print_step == 0: print(n_mol_processed, "mols processed in block",file_num,"of size", n_mols)
        # for i, query_fp in enumerate(query_fp_matrix):
        for i in numba.prange(query_fp_matrix.shape[0]):
            query_fp = query_fp_matrix[i]
            # print(numba.typeof(query_fp))
            # print(numba.typeof(finPrint))
            simil = compute_2FingerPrints_similarity(query_fp, bin_finPrint)
            # print(compute_2FingerPrints_similarity.signatures)
            # print(np.sum(query_fp), np.sum(finPrint), simil)
            less_similar_idx = np.argmin(matched_similarities[i, :])

            if simil > matched_similarities[i, less_similar_idx]:
                matched_similarities[i, less_similar_idx] = simil
                # matched_ids[i, less_similar_idx, :] = [file_num, j//fp_bipts]
                matched_ids[i, less_similar_idx, 0] = file_num
                matched_ids[i, less_similar_idx, 1] = j//fp_bipts


    return matched_similarities, matched_ids


def process_one_subFile_numba(query_fps_mat, fileNum_chunkFname, n_hits_per_smi):

    file_num, chunk_fname = fileNum_chunkFname
    with open(chunk_fname, "rb") as f:
        db_many_fps_bool = decompressFingerprint_npStr (f.read() )

    matched_similarities, matched_ids = compute_FingerPrints_similarity_from_allbools(query_fps_mat, db_many_fps_bool, file_num, n_hits_per_smi)

    return  matched_similarities, matched_ids


def process_one_subFile_numpy(query_fps_mat, fileNum_chunkFname, n_hits_per_smi): #this is slower than numba for small number of queries but A LOT faster large query numbers . Consumes a lot of memory # TODO use dask array to use mapped arrays

    #TODO: Add gpu support with https://docs.cupy.dev/en/stable/tutorial/basic.html
    #TODO: estimante memory requirements and process chunk by chunk

    file_num, chunk_fname = fileNum_chunkFname

    with open(chunk_fname, "rb") as f:
        db_many_fps_bool = decompressFingerprint_npStr (f.read() )

    sim_matrix = jaccard_vectorized(query_fps_mat, np.frombuffer(db_many_fps_bool, dtype=np.bool).reshape(-1, query_fps_mat.shape[1]))
    max_sim_idx_per_compound = np.argsort(-sim_matrix, axis=1)[:, :n_hits_per_smi]
    matched_similarities = np.ones((len(query_fps_mat), n_hits_per_smi)) * -1
    matched_similarities[:,:max_sim_idx_per_compound.shape[1]] = np.take_along_axis(sim_matrix, max_sim_idx_per_compound, axis=1)
    matched_ids = -1*np.ones( matched_similarities.shape+(2,) )
    matched_ids[...,0] *= -file_num
    matched_ids[:,:max_sim_idx_per_compound.shape[1], 1] = max_sim_idx_per_compound
    return  matched_similarities, matched_ids


def combine_two_chunk_searchs(cs1, cs2):
    '''

    :param cs1: tuple( similarities, ids)
    :param cs2: same as cs1
    :return:
    '''

    out_simil, out_ids = cs1[0].copy(), cs1[1].copy()
    biggerSim_second_mask = cs1[0] < cs2[0]

    out_simil[biggerSim_second_mask] = cs2[0][biggerSim_second_mask]
    out_ids[biggerSim_second_mask] = cs2[1][biggerSim_second_mask]
    return out_simil, out_ids

def combine_search_jsons(path_to_jsons):

    def combine_two_dict(d1, d2, ):
        new_d = {}
        assert len(d1)==len(d2), "Error, both dicts should be equaly-sized"
        for query_smi, matches_list1 in d1.items():
            matches_list2 = d2[query_smi]
            assert len(matches_list1) == len(matches_list2), "Error, both list should be equaly-sized"
            final_matches = []
            for match1, match2 in zip(matches_list1, matches_list2):
                similarity1 = match1[0]
                similarity2 = match2[0]
                if similarity1 > similarity2:
                    final_matches.append( match1)
                else:
                    final_matches.append( match2)
            new_d[query_smi] = final_matches
        return new_d

    def load_json(fname):
        with open(fname) as f:
            return json.load(f)

    filenames = [os.path.join(path_to_jsons, fname) for fname in os.listdir(path_to_jsons) ]
    bag = db.from_sequence( filenames[1:]).map(load_json).fold(combine_two_dict, initial= load_json(filenames[0]) )
    final_search= bag.compute()
    print(final_search)
    return final_search

def search_smi_list(query_smi_list, database_dir, n_hits_per_smi=30, output_name=None, verbose=True):

    starting_time = time.time()

    def compute_query_fp( qsmi):
        fp = get_fingerPrint_as_npBool(qsmi)
        if fp is None:
            return None
        else:
            return qsmi, fp

    query_fps = db.from_sequence(query_smi_list).map(compute_query_fp).filter(None.__ne__)

    query_smi_list, query_fps,  = zip(*query_fps.compute())
    query_smi_list = list(query_smi_list)
    query_fps = list(query_fps)

    query_fps = np.stack(query_fps, axis=0) # N_queriesxfingerprint_n_bits type boolean (wasting 7/8 memory but faster)



    matched_similarities= np.ones( (len(query_fps), n_hits_per_smi) ) * -1              #query_id, hit_num, similarity
    matched_ids = np.ones( (len(query_fps), n_hits_per_smi, 2), dtype= np.int64 ) * -1  #query_id, hit_num, [ file_id, hit_id]

    if query_fps.shape[0] > 5:
        process_one_subFile = lambda fname: process_one_subFile_numpy(query_fps, fname, n_hits_per_smi)
    else:
        process_one_subFile = lambda fname: process_one_subFile_numba(query_fps, fname, n_hits_per_smi) #faster for small number of query mols

    fingerprints_dir = os.path.join(database_dir, "fingerprints")
    filenames = filter( lambda x:  x.endswith(".fingerprints.BitVect"), sorted(os.listdir(fingerprints_dir)))
    filenames =  list( map( lambda x:  os.path.join(fingerprints_dir,x), filenames ) )


    bag = db.from_sequence( enumerate(filenames)).map(process_one_subFile).fold(combine_two_chunk_searchs, initial=(matched_similarities, matched_ids))
    matched_similarities, matched_ids = bag.compute()

    # print( matched_similarities )
    # print( matched_ids )

    compounds_name = os.path.join(database_dir, "compounds.sqlite")


    basenames = [ os.path.basename(fname).split(".")[0] for fname in filenames]
    # matches_db_entries_compound = [ (int(rowNum), basenames[fileNum]) for fileNum, rowNum in matched_ids.reshape(-1, 2) ]

    con = sqlite3.connect(compounds_name)
    cur = con.cursor()

    resultsDict = OrderedDict([])
    n_found = 0
    for queryIdx in range(matched_ids.shape[0]):
        matches_db_entries_compound = [(int(rowNum), basenames[fileNum]) for fileNum, rowNum in matched_ids[queryIdx,...].reshape(-1, 2) if rowNum>=0]
        matches_sim_compound = matched_similarities[queryIdx,...]
        query_smi = query_smi_list[queryIdx]
        if verbose: print("> Query: %s"%query_smi)
        # print( matches_db_entries_compound )
        if len(matches_db_entries_compound) > 0:
            resultsDict[query_smi] = []
            found_in_db=0
            for sim, match_entry in sorted(zip(matches_sim_compound, matches_db_entries_compound), key= lambda x: -x[0]):
                # print( match_entry )
                for row in cur.execute("SELECT compoundId, smi FROM smiles NATURAL JOIN ("
                                       " SELECT compoundId FROM compounds WHERE rowNum=? AND fileSource=?)", match_entry):
                    resultsDict[query_smi].append((sim, row[0], row[1]))
                    if verbose: print("%.4f\t%s\t%s" % (sim, row[0], row[1]))
                    n_found +=1
                    found_in_db+=1
            assert  found_in_db == len(matches_db_entries_compound), "Error in database!!"
        else:
            if verbose: print("No matching compounds were found")

    con.close()
    total_time = time.time() - starting_time
    if verbose: print("Total time for %d smi: %s" % ( n_found, str(datetime.timedelta(seconds=total_time))))
    for key in resultsDict:
        resultsDict[key].sort(key=lambda x: -x[0])

    if output_name:
        with open(output_name, "w") as f:
            json.dump(resultsDict, f)  # >>> data = json.loads(json_str, object_pairs_hook=OrderedDict)
    return resultsDict


def mainSearch():
    import sys
    import argparse
    parser = ArgumentParser(prog="fast_similarity_search", description="Find the K most similar compounds in the database")
    parser.add_argument('smiles_query', nargs=None, type=argparse.FileType('r'), default=sys.stdin,
                        help="smiles file with as many smiles as rows")

    parser.add_argument('-d', '--database_dir', help="the directory where compounds database was compiled",  required=True)

    parser.add_argument('-m', '--metric', nargs=None, choices=["Tanimoto"], default="Tanimoto", required=False,
                        help="metric to use")

    parser.add_argument('-n', '--n_hits_per_smi', type=int, default=30,
                        help="K highest score molecules to retrieve per query smiles ")

    parser.add_argument('-o', '--output_name', type=str, required=True,
                        help="The fname for a json file where search results will be stored")

    parser.add_argument('-v', '--verbose', action="store_true", default=False,
                        help="Print to stdout working information ")

    args = parser.parse_args()
    query_smi_list =  args.smiles_query.read().splitlines()
    query_smi_list = [ x.strip() for x in query_smi_list]

    assert len(query_smi_list) >0, "Error, no smiles provided"

    dask_client = get_parallel_client()

    # numba.set_num_threads(multiprocessing.cpu_count()/n_dask_workers)
    search_smi_list(query_smi_list, args.database_dir, n_hits_per_smi=args.n_hits_per_smi, output_name=args.output_name,
                    verbose=args.verbose)

    dask_client.shutdown()

def testSearch():
    fp_outdir = "/home/ruben/oxford/enamine/fingerprints"
    save_fname = "/home/ruben/tmp/mols/first_partition.json"
    query_smi_list = ["CC1CCSCCN1CCC1=CN=CC(F)=C1", "COC", "C1C(C(C(C(O1)O)O)O)O"]
    found = search_smi_list(query_smi_list, fp_outdir, n_hits_per_smi=3, output_name= save_fname)
    print( found )

def testCombine():
    jsons_dir = "/home/ruben/tmp/mols/"
    res = combine_search_jsons(jsons_dir)

if __name__ == "__main__":
    mainSearch()

    '''
echo -e "CC1CCSCCN1CCC1=CN=CC(F)=C1\nCOC\nC1C(C(C(C(O1)O)O)O)O" | N_CPUS=2 python -m fragmenstein.external.enamine_realDB_search.similarity_fast_searcher -d /home/ruben/oxford/enamine/fingerprints -o /home/ruben/tmp/mols/first_partition.json -

tail  -n 1 ~/oxford/enamine/cxsmiles/big_example1.txt | cut -f 1 | N_CPUS=1 python -m fragmenstein.external.enamine_realDB_search.similarity_fast_searcher -d /home/ruben/oxford/enamine/fingerprints -o /home/ruben/tmp/mols/first_partition.json -v -


    '''