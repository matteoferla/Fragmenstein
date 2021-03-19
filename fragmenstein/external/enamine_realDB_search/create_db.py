import datetime
import os
import sqlite3
import subprocess
import tempfile
import numpy as np
import pandas as pd
import threading
import time
from queue import Queue

import dask.bag as db
import dask.dataframe as dd
from dask.distributed import futures_of, as_completed
from joblib import Parallel, delayed

from fragmenstein.external.enamine_realDB_search.common import  computeFingerprint_np_Str
from fragmenstein.utils.config_manager import ConfigManager
from fragmenstein.utils.parallel_utils import get_parallel_client



N_LINES_PER_CHUNK = int(5e5) # int(1e6)



def chunk_fname(fname,  chunked_dir, n_lines_per_chunk=N_LINES_PER_CHUNK):

    fname_base = os.path.basename(fname).split(".")[0]
    chunk_prefix = os.path.join(chunked_dir, fname_base)

    decompression = ""
    if fname.endswith(".bz2"):
        decompression =  "  | lbzip2 -dc " # lbzip2 is way faster than regular bzip2
        if ConfigManager.N_CPUS > 1:
            decompression += "-n% d "%ConfigManager.N_CPUS

    cmd = "cat "+fname+ decompression+" | tail -n +2  | split -l "+str(n_lines_per_chunk)+" - "+ chunk_prefix # +"--filter='lbzip2 > $FILE.csv.gz' "
    # input( cmd)
    proc = subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None,
                                  close_fds=True)

    proc.wait() #todo, check returncode

    return list(filter(os.path.isfile, (os.path.join(chunked_dir, name) for name in os.listdir(chunked_dir)) ) )

    #SPOOLING DOES NOT MAKE SENSE IF USING dask bag, as it converts iterable to list
    #SPOOLING_TIME =5
    # keep_spooling = True
    # available_paths = Queue()
    # already_seen_names = set([])
    #
    # def find_unseen_names():
    #     all_fnames = filter(os.path.isfile, (os.path.join(chunked_dir, name) for name in os.listdir(chunked_dir)))
    #     # remove already seen
    #     all_fnames = set(all_fnames) - already_seen_names
    #     unseen_names = sorted(all_fnames, key=os.path.getmtime)
    #     return unseen_names
    #
    # def spool():
    #     while keep_spooling:
    #         paths = find_unseen_names()
    #         if len(paths)>1:
    #             for path in paths[:-1]:
    #                 available_paths.put( path )
    #                 already_seen_names.add( path )
    #         time.sleep(SPOOLING_TIME)
    #
    # t1 = threading.Thread(target=spool, daemon=True)
    # t1.start()
    #
    #
    # while proc.poll() is None: #while still compressing
    #     while not available_paths.empty():
    #         fname = available_paths.get()
    #         yield fname
    #     time.sleep(SPOOLING_TIME+1)
    # else: #afther while loop, will the the thread
    #     keep_spooling = False
    #
    # t1.join()
    #
    # while not available_paths.empty():
    #     fname = available_paths.get()
    #     yield fname
    #
    # for fname in find_unseen_names():
    #     yield fname


def process_cxsmi_file( fname, compunds_db_fname, binaries_dir, n_lines_per_chunk=N_LINES_PER_CHUNK):

    starting_time = time.time()

    print( fname )
    with tempfile.TemporaryDirectory() as tmpdir:
        print(tmpdir)

        chunked_names = chunk_fname(fname, tmpdir, n_lines_per_chunk=n_lines_per_chunk)

        def process_one_chunkedFile(chunked_fname):

            chunked_basename = os.path.basename(chunked_fname).split(".")[0]

            print("processing %s "%chunked_basename, flush=True)

            data = dd.read_csv(chunked_fname, sep="\t", header=None, usecols=[0,1])
            print("%s raw data read"%chunked_basename, flush=True)

            t = time.time()

            binary_name = os.path.join(binaries_dir, chunked_basename+".fingerprints.BitVect")

            data = data.compute()

            ids = []
            smis = []
            n_fingerprints =0
            with open(binary_name, "wb") as bin_f:
                for row in data.itertuples():
                    __, smi, cid = row[:3]
                    fp = computeFingerprint_np_Str(smi)
                    if fp is not None:
                        bin_f.write(fp)
                        smis.append( smi)
                        ids.append(cid )
                        n_fingerprints+=1

            id_table=  pd.DataFrame( (cid, chunked_basename, i) for i,cid in enumerate(ids) )
            id_table.columns = ["compoundId", "fileSource", "rowNum"]

            smis_table = pd.DataFrame(zip(ids, smis))
            smis_table.columns = ["compoundId", "smi" ]

            print("%s fingenprints computed (%s s)" % (chunked_basename, time.time() - t), flush=True)

            return chunked_basename, id_table, smis_table

        b = db.from_sequence(chunked_names).map( process_one_chunkedFile )
        b = b.persist()
        futures_b = futures_of(b)

        con = sqlite3.connect( compunds_db_fname)
        total_count = 0
        for fut in as_completed(futures_b):
            for res in fut.result():
                (chunked_basename, compounds_df, smiles_df) = res
                print( compounds_df.query("rowNum==0"))
                compounds_df.to_sql("compounds", con, index=False, if_exists="append")
                smiles_df.to_sql("smiles", con, index=False, if_exists="append")
                con.commit()
                total_count+= len(compounds_df)
                print("%s was commited. Processed smis: %d"%(chunked_basename, total_count), flush=True)

        con.commit()
        con.close()
    total_time = time.time() - starting_time
    print( "Total time for %s ( %d smi): %s"%(fname, total_count, str(datetime.timedelta(seconds=total_time)) ))

def process_all_files(cxsmiles_dir, outdir, n_lines_per_chunk=N_LINES_PER_CHUNK, work_in_memory=False, *args, **kwargs):

    actual_compunds_db_fname = os.path.join( outdir, "compounds.sqlite")
    assert  not os.path.exists(actual_compunds_db_fname), "Error, sqlite file %s already existing"%actual_compunds_db_fname

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    with tempfile.TemporaryDirectory(dir="/dev/shm", suffix="_"+str(abs(hash(cxsmiles_dir)%2048))) as tmp:
        if work_in_memory:
            compunds_db_fname = os.path.join(tmp, os.path.basename(actual_compunds_db_fname) )
        else:
            compunds_db_fname = actual_compunds_db_fname


        con = sqlite3.connect( compunds_db_fname)

        cur = con.cursor()

        cur.execute('''CREATE TABLE compounds
                       ( compoundId VARCHAR(20) PRIMARY KEY, fileSource VARCHAR(50), rowNum INTEGER)''')

        cur.execute('''CREATE INDEX index_on_file_row ON compounds(fileSource, rowNum)
        ''')

        cur.execute('''CREATE TABLE smiles
                       (compoundId VARCHAR(20) PRIMARY KEY, smi TEXT)''')
        con.commit()


        if os.path.isdir(cxsmiles_dir):
            fnames = map(lambda x: os.path.join(cxsmiles_dir, x), os.listdir(cxsmiles_dir))
        else:
            fnames = [cxsmiles_dir]

        binaries_dir = os.path.join( outdir, "fingerprints")
        if not os.path.exists( binaries_dir ):
            os.mkdir( binaries_dir )

        Parallel(n_jobs=1)(delayed(process_cxsmi_file)(fname,compunds_db_fname, binaries_dir,
                                                       n_lines_per_chunk)  for fname in fnames)


        if work_in_memory:
            print( "Dumping db to disk")

            bck = sqlite3.connect(actual_compunds_db_fname)
            with bck:
                con.backup(bck, pages=-1)
            bck.close()

        con.close()

# def main():
#     import sys
#     cxsmiles_dir = sys.argv[1]
#     fp_outdir = sys.argv[2]
#
#     dask_client = get_parallel_client()
#     process_all_files(cxsmiles_dir, fp_outdir)
#     dask_client.shutdown()


def main():
    import sys
    import argparse
    from fragmenstein.utils.cmd_parser import ArgumentParser
    parser = ArgumentParser(prog="compile_similarity_search_db", description="Compiles a database for similarity search on cxsmiles files coming from enamine")


    parser.add_argument('-i', '--cxsmiles_dir', help="path to one single cssmiles file or a dicectory containing several",  required=True)

    parser.add_argument('-o', '--outdir', type=str, required=True,
                        help="The fname for a json file where search results will be stored")

    parser.add_argument('-f', '--fingerprint', nargs=None, choices=["Morgan"], default="Morgan", required=False,
                        help="fingerprint to use")

    parser.add_argument( '--work_in_memory', action="store_true", default=False,
                        help="if work in memory, database will be created in memory and dumped to file after completion")

    parser.add_argument('-n', '--n_lines_to_chunk', type=int, default=N_LINES_PER_CHUNK,
                        help="Number of lines to chunkenize input files for parallel processing. Default: %(default)s")


    # parser.add_argument('-v', '--verbose', action="store_true", default=False,
    #                     help="Print to stdout working information ")

    args = parser.parse_args()
    print(args)
    dask_client = get_parallel_client()
    process_all_files(**vars(args))
    dask_client.shutdown()

def testCreate():
    cxsmiles_dir = "/home/ruben/oxford/enamine/cxsmiles"  # Enamine_REAL_HAC_21_22_CXSMILES.cxsmiles.bz2"
    fp_outdir = "/home/ruben/oxford/enamine/fingerprints"
    process_all_files(cxsmiles_dir, fp_outdir)

if __name__ == "__main__":
    # testCreate()
    main()


    '''

N_CPUS=4 python -m fragmenstein.external.enamine_realDB_search.create_db -i /home/ruben/oxford/enamine/cxsmiles -o /home/ruben/oxford/enamine/fingerprints

python -m fragmenstein.external.condor_queue.send_to_condor --env_vars EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json PATH=/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin:$PATH DASK_WORKER_MEMORY=4GB N_CPUS=42 --ncpus 42 "/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python -m fragmenstein.external.enamine_realDB_search.create_db  --work_in_memory -i /data/xchem-fragalysis/sanchezg/oxford/enamine/full_cxsmiles/full_cxsmiles/Enamine_REAL_HAC_21_22_CXSMILES.cxsmiles.bz2 -o /data/xchem-fragalysis/sanchezg/oxford/enamine/fingerprints_db/Enamine_REAL_HAC_21_22_CXSMILES"


    '''