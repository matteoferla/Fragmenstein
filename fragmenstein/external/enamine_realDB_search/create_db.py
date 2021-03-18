import datetime
import os
import sqlite3
import subprocess
import tempfile
import threading
import time
from queue import Queue

import dask.bag as db
import dask.dataframe as dd
from dask.distributed import futures_of, as_completed
from joblib import Parallel, delayed

from fragmenstein.external.enamine_realDB_search.common import computeFingerprintStr, computeFingerprint_np_Str
from fragmenstein.utils.config_manager import ConfigManager
from fragmenstein.utils.parallel_utils import get_parallel_client



N_LINES_PER_CHUNK = int(5e5) # int(1e6)

SPOOLING_TIME = 5


def chunk_fname(fname,  chunked_dir):

    fname_base = os.path.basename(fname).split(".")[0]
    chunk_prefix = os.path.join(chunked_dir, fname_base)

    decompression = ""
    if fname.endswith(".bz2"):
        decompression =  "  | lbzip2 -dc " # lbzip2 is way faster than regular bzip2
        if ConfigManager.N_CPUS > 1:
            decompression += "-n% d "%ConfigManager.N_CPUS

    cmd = "tail -n +2  "+fname+decompression+" | split -l "+str(N_LINES_PER_CHUNK)+" - "+ chunk_prefix # +"--filter='lbzip2 > $FILE.csv.gz' "
    # input( cmd)
    proc = subprocess.Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None,
                                  close_fds=True)
    keep_spooling = True
    available_paths = Queue()
    already_seen_names = set([])

    def find_unseen_names():
        all_fnames = filter(os.path.isfile, (os.path.join(chunked_dir, name) for name in os.listdir(chunked_dir)))
        # remove already seen
        all_fnames = set(all_fnames) - already_seen_names
        unseen_names = sorted(all_fnames, key=os.path.getmtime)
        return unseen_names

    def spool():
        while keep_spooling:
            paths = find_unseen_names()
            if len(paths)>1:
                for path in paths[:-1]:
                    available_paths.put( path )
                    already_seen_names.add( path )
            time.sleep(SPOOLING_TIME)

    t1 = threading.Thread(target=spool, daemon=True)
    t1.start()


    while proc.poll() is None: #while still compressing
        while not available_paths.empty():
            fname = available_paths.get()
            yield fname
        time.sleep(SPOOLING_TIME+1)
    else: #afther while loop, will the the thread
        keep_spooling = False

    t1.join()

    while not available_paths.empty():
        fname = available_paths.get()
        yield fname

    for fname in find_unseen_names():
        yield fname


def process_cxsmi_file( fname,  outdir):

    starting_time = time.time()

    print( fname )
    with tempfile.TemporaryDirectory() as tmpdir:
        print(tmpdir)

        chunked_names = chunk_fname(fname, tmpdir)

        compounds_name = os.path.join( outdir, "compounds.sqlite")
        binaries_dir = os.path.join( outdir, "fingerprints")
        if not os.path.exists( binaries_dir ):
            os.mkdir( binaries_dir )

        def process_one_chunkedFile(chunked_fname):

            chunked_basename = os.path.basename(chunked_fname).split(".")[0]

            print("processing %s "%chunked_basename, flush=True)

            data = dd.read_csv(chunked_fname, sep="\t", header=None, usecols=[0,1])
            print("%s raw data read"%chunked_basename, flush=True)

            t = time.time()

            binary_name = os.path.join(binaries_dir, chunked_basename+".fingerprints.BitVect")

            with open(binary_name, "wb") as bin_f:
                for smi in data[0]:
                    fp = computeFingerprint_np_Str(smi)
                    if fp is not None:
                        bin_f.write(fp)

            print(time.time() - t)
            print("%s fingenprints computed" % chunked_basename, flush=True)

            data = data.compute()
            id_table=  data[[1]]
            id_table.columns = ["compoundId"]
            id_table["fileSource"] = chunked_basename
            id_table["rowNum"] = id_table.index
            smis_table = data[[1,0]]
            smis_table.columns = ["compoundId", "smi" ]
            print("%s processed"%chunked_basename, flush=True)

            return id_table, smis_table

        b = db.from_sequence(chunked_names).map( process_one_chunkedFile )
        b = b.persist()
        futures_b = futures_of(b)

        con = sqlite3.connect( compounds_name)
        total_count = 0
        for fut in as_completed(futures_b):
            for res in fut.result():
                (compounds_df, smiles_df) = res
                compounds_df.to_sql("compounds", con, index=False, if_exists="append")
                smiles_df.to_sql("smiles", con, index=False, if_exists="append")
                con.commit()
                total_count+= len(compounds_df)
                print("processed smis: %d"%total_count)

        con.commit()
        con.close()
    total_time = time.time() - starting_time
    print( "Total time for %s ( %d smi): %s"%(fname, total_count, str(datetime.timedelta(seconds=total_time)) ))

def process_all_files(inpath, outdir):

    compounds_name = os.path.join( outdir, "compounds.sqlite")
    assert  not os.path.exists(compounds_name), "Error, sqlite file %s already existing"%compounds_name
    con = sqlite3.connect( compounds_name)
    cur = con.cursor()

    cur.execute('''CREATE TABLE compounds
                   ( compoundId VARCHAR(20) PRIMARY KEY, fileSource VARCHAR(50), rowNum INTEGER)''')

    cur.execute('''CREATE INDEX index_on_file_row ON compounds(fileSource, rowNum)
    ''')

    cur.execute('''CREATE TABLE smiles
                   (compoundId VARCHAR(20) PRIMARY KEY, smi TEXT)''')
    con.commit()

    if os.path.isdir(inpath):
        fnames = map( lambda x: os.path.join(inpath, x), os.listdir(inpath))
    else:
        fnames = [inpath]
    Parallel(n_jobs=1)(delayed(process_cxsmi_file)(fname, outdir)  for fname in fnames )



def merge_squlites():
    import sqlite3
    con3 = sqlite3.connect("combine.db")

    con3.execute("ATTACH 'results_a.db' as dba")

    con3.execute("BEGIN")
    for row in con3.execute("SELECT * FROM dba.sqlite_master WHERE type='table'"):
        combine = "INSERT INTO " + row[1] + " SELECT * FROM dba." + row[1]
        print(combine)
        con3.execute(combine)
    con3.commit()
    con3.execute("detach database dba")
    raise NotImplementedError("This is only an example")

def testCreate():
    cxsmiles_dir = "/home/ruben/oxford/enamine/cxsmiles" #Enamine_REAL_HAC_21_22_CXSMILES.cxsmiles.bz2"
    fp_outdir = "/home/ruben/oxford/enamine/fingerprints"
    process_all_files(cxsmiles_dir, fp_outdir)

def main():
    import sys
    cxsmiles_dir = sys.argv[1]
    fp_outdir = sys.argv[2]

    dask_client = get_parallel_client()
    process_all_files(cxsmiles_dir, fp_outdir)
    dask_client.shutdown()

if __name__ == "__main__":
    # testCreate()
    main()


    '''

N_CPUS=4 python -m fragmenstein.external.enamine_realDB_search.create_db /home/ruben/oxford/enamine/cxsmiles /home/ruben/oxford/enamine/fingerprints

    '''