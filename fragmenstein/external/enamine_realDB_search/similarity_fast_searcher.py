import datetime
import json
import os
import sqlite3
import time
from collections import OrderedDict

import dask.bag as db
import numpy as np
from rdkit import DataStructs
from rdkit.DataStructs import cDataStructs

from fragmenstein.external.enamine_realDB_search.common import FINGERPRINT_NBITS, get_fingerPrint
from fragmenstein.utils.parallel_utils import get_parallel_client


def compute_2FingerPrints_similarity(fp1, fp2):
    return DataStructs.FingerprintSimilarity(fp1, fp2)


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

def search_smi_list(query_smi_list, db_files_dir, n_hits_per_smi=30, save_fname=None, verbose=True):

    starting_time = time.time()

    query_fps = db.from_sequence(query_smi_list).map(get_fingerPrint) #.filter(None.__ne__)
    query_fps = query_fps.compute()

    chunk_bytes = FINGERPRINT_NBITS // 8

    def process_one_chunk(fileNum_chunkFname):

        file_num, chunk_fname = fileNum_chunkFname

        matched_similarities = np.ones((len(query_fps), n_hits_per_smi)) * -1  # query_id, hit_num, similarity
        matched_ids = np.ones((len(query_fps), n_hits_per_smi, 2),
                              dtype=np.int64) * -1  # query_id, hit_num, [ file_id, hit_id]

        with open(chunk_fname, "rb") as f:
            bin_finPrint = f.read(chunk_bytes)
            chunk_number = 0
            while bin_finPrint:
                # process
                finPrint = cDataStructs.CreateFromBinaryText(bin_finPrint)
                for i, query_fp in enumerate(query_fps):
                    simil = compute_2FingerPrints_similarity(query_fp, finPrint)
                    less_similar_idx = np.argsort(matched_similarities[i, :])[0]
                    if simil > matched_similarities[i, less_similar_idx]:
                        matched_similarities[i, less_similar_idx] = simil
                        matched_ids[i, less_similar_idx, :] = [file_num, chunk_number]

                chunk_number += 1
                bin_finPrint = f.read(chunk_bytes)

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

    matched_similarities= np.ones( (len(query_fps), n_hits_per_smi) ) * -1              #query_id, hit_num, similarity
    matched_ids = np.ones( (len(query_fps), n_hits_per_smi, 2), dtype= np.int64 ) * -1  #query_id, hit_num, [ file_id, hit_id]

    fingerprints_dir = os.path.join(db_files_dir, "fingerprints")
    filenames = filter( lambda x:  x.endswith(".fingerprints.BitVect"), sorted(os.listdir(fingerprints_dir)))
    filenames =  list( map( lambda x:  os.path.join(fingerprints_dir,x), filenames ) )


    bag = db.from_sequence( enumerate(filenames)).map(process_one_chunk).fold(combine_two_chunk_searchs, initial=(matched_similarities, matched_ids))
    matched_similarities, matched_ids = bag.compute()

    # print( matched_similarities )
    # print( matched_ids )


    compounds_name = os.path.join( db_files_dir, "compounds.sqlite")


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
        print("> Query: %s"%query_smi)
        print( matches_db_entries_compound )
        if len(matches_db_entries_compound) > 0:
            resultsDict[query_smi] = []
            for sim, match_entry in zip(matches_sim_compound, matches_db_entries_compound):
                for row in cur.execute("SELECT compoundId, smi FROM smiles NATURAL JOIN ("
                                       " SELECT compoundId FROM compounds WHERE rowNum=? AND fileSource=?)", match_entry):
                    resultsDict[query_smi].append((sim, row[0], row[1]))
                    if verbose: print("%.4f\t%s\t%s" % (sim, row[0], row[1]))
                    n_found +=1
        else:
            if verbose: print("No matching compounds were found")

    con.close()
    total_time = time.time() - starting_time
    if verbose: print("Total time for %d smi: %s" % ( n_found, str(datetime.timedelta(seconds=total_time))))
    for key in resultsDict:
        resultsDict[key].sort(key=lambda x: -x[0])

    if save_fname:
        with open(save_fname, "w") as f:
            json.dump(resultsDict, f)  # >>> data = json.loads(json_str, object_pairs_hook=OrderedDict)
    return resultsDict

def testSearch():
    fp_outdir = "/home/ruben/oxford/enamine/fingerprints"
    save_fname = "/home/ruben/tmp/mols/first_partition.json"
    query_smi_list = ["CC1CCSCCN1CCC1=CN=CC(F)=C1", "COC", "C1C(C(C(C(O1)O)O)O)O"]
    found = search_smi_list(query_smi_list, fp_outdir, n_hits_per_smi=3, save_fname= save_fname)
    print( found )

def testCombine():
    jsons_dir = "/home/ruben/tmp/mols/"
    res = combine_search_jsons(jsons_dir)

if __name__ == "__main__":
    dask_client = get_parallel_client()
    # testSearch()
    testCombine()
    dask_client.shutdown()

    '''
python -m fragmenstein.external.enamine_realDB_search.similarity_fast_searcher

    '''