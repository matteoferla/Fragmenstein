import os
import sqlite3
import time, datetime
from itertools import chain

import dask.bag as db
import numpy as np
from joblib import Parallel, delayed, load
from rdkit import DataStructs
from rdkit.DataManip.Metric import rdMetricMatrixCalc
from rdkit.DataStructs import cDataStructs

from fragmenstein.external.enamine_realDB_search.fingerprint_params import FINGERPRINT_NBITS, get_fingerPrint
from fragmenstein.utils.parallel_utils import get_parallel_client


def compute_2FingerPrints_similarity(fp1, fp2):
    return DataStructs.FingerprintSimilarity(fp1, fp2)


# def search_smi_list_oneFilel(query_fps, db_oneFile):
#
#     db_oneFile_fps = load( db_oneFile)
#     def compute_similarity(fp):
#         dist_mat = rdMetricMatrixCalc.GetTanimotoDistMat( chain.from_iterable([db_oneFile_fps, [fp]]) )
#
#     dists = Parallel(n_jobs=1)(delayed(get_fingerPrint)(fp) for fp in query_fps)


def search_smi_list_global(query_smi_list, db_files_dir,  n_hits_per_smi=30):

    starting_time = time.time()

    query_fps = db.from_sequence(query_smi_list).map(get_fingerPrint) #.filter(None.__ne__)
    query_fps = query_fps.compute()

    matched_similarities= np.ones( (len(query_fps), n_hits_per_smi) ) * -1              #query_id, hit_num, similarity
    matched_ids = np.ones( (len(query_fps), n_hits_per_smi, 2), dtype= np.int64 ) * -1  #query_id, hit_num, [ file_id, hit_id]

    fingerprints_dir = os.path.join(db_files_dir, "fingerprints")
    filenames = list( filter( lambda x: x.endswith(".fingerprints.BitVect"), sorted(os.listdir(fingerprints_dir)) ))
    chunk_bytes = FINGERPRINT_NBITS // 8
    for file_num, fname in enumerate(filenames): #TODO. Ensure fingerprints files always preserve order
        with open( os.path.join(fingerprints_dir, fname), "rb") as f:
            bin_finPrint = f.read(chunk_bytes)
            chunk_number = 0
            while bin_finPrint:
                #process
                finPrint = cDataStructs.CreateFromBinaryText(bin_finPrint)
                for i, query_fp in enumerate(query_fps):
                    simil = compute_2FingerPrints_similarity(query_fp, finPrint)
                    less_similar_idx = np.argsort(matched_similarities[i,:])[0]
                    if simil > matched_similarities[i,less_similar_idx]:
                        matched_similarities[i, less_similar_idx] = simil
                        matched_ids[i, less_similar_idx, :] = [ file_num, chunk_number]

                chunk_number +=1
                bin_finPrint = f.read(chunk_bytes)


    print( matched_similarities )
    print( matched_ids )


    compounds_name = os.path.join( db_files_dir, "compounds.sqlite")


    basenames = [ fname.split(".")[0] for fname in filenames]
    # matches_db_entries_compound = [ (int(rowNum), basenames[fileNum]) for fileNum, rowNum in matched_ids.reshape(-1, 2) ]

    con = sqlite3.connect(compounds_name)
    cur = con.cursor()

    for queryIdx in range(matched_ids.shape[0]):
        matches_db_entries_compound = [(int(rowNum), basenames[fileNum]) for fileNum, rowNum in matched_ids[queryIdx,...].reshape(-1, 2) if rowNum>=0]
        matches_sim_compound = matched_similarities[queryIdx,...]
        print("> Query: %s"%query_smi_list[queryIdx])
        print( matches_db_entries_compound )
        if len(matches_db_entries_compound) > 0:

            for sim, match_entry in zip(matches_sim_compound, matches_db_entries_compound):
                for row in cur.execute("SELECT compoundId, smi FROM smiles NATURAL JOIN ("
                                       " SELECT compoundId FROM compounds WHERE rowNum=? AND fileSource=?)", match_entry):
                    print("%s\t%.4f\t%s"%(row[0], sim, row[1]) )
        else:
            print("No matching compounds were found")

    con.close()
    total_time = time.time() - starting_time
    print( "Total time for %s ( %d smi): %s"%(fname, i, str(datetime.timedelta(seconds=total_time)) ))

        # for file_num, fname in enumerate(filenames):
    #     matches_in_file_mask = matched_ids[..., 0] == file_num
    #     matching_compound_lines = matched_ids[matches_in_file_mask, ..., 1]
    #     matching_similarites = matched_similarities[matches_in_file_mask, ...]
    #
    #     print(matching_compound_lines)
    #     print(matching_similarites)
    #
    #     if len(matching_similarites)>0:
    #
    #         cur.executemany("SELECT compoundId FROM compounds WHERE fileSource=? and rowNum=?", )
    #         compounds_index_fname = os.path.join(db_files_dir, fname.replace(".fingerprints.BitVect", ".compounds.tab"))
    #         with open(compounds_index_fname) as f: #TODO: write it as a binary file of fixed size to use seek
    #             for i, line in enumerate(f):
    #                 if i in matching_compound_lines:
    #                     print(line )
    # input(query_fps)



def testSearch():
    fp_outdir = "/home/ruben/oxford/enamine/fingerprints"
    query_smi_list = ["CC1CCSCCN1CCC1=CN=CC(F)=C1", "COC", "C1C(C(C(C(O1)O)O)O)O"]
    search_smi_list_global(query_smi_list, fp_outdir, n_hits_per_smi=3)

if __name__ == "__main__":
    dask_client = get_parallel_client()
    testSearch()

    '''
python -m fragmenstein.external.enamine_realDB_search.similarity_fast_searcher

    '''