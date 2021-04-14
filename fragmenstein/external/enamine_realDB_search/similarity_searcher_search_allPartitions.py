import os
import tempfile
import time

from subprocess import check_call

import re

WAIT_TIME_LAUNCH_QUEUE = 5
SHARED_TMP = "/data/xchem-fragalysis/sanchezg/tmp/"

def findDB_partitions(db_path):
  partitions = []
  for fname in os.listdir(db_path):
    full_path = os.path.join(db_path, fname)
    if os.path.isfile(os.path.join(full_path, "compounds.sqlite")):
      partitions.append( full_path )

  return partitions



def parse_memsize(size):
  #https://stackoverflow.com/questions/42865724/parse-human-readable-filesizes-into-bytes/42865957#42865957
  units = {"KB": 2 ** -10, "MB": 1, "GB": 2 ** 10, "TB": 2 ** 20} #UNIT_TO_MBs
  size = size.upper()
  #print("parsing size ", size)
  if not re.match(r' ', size):
      size = re.sub(r'([KMGT]?B)', r' \1', size)
  number, unit = [string.strip() for string in size.split()]
  return int(float(number)*units[unit])

def launch_searcher(run_locally=False, **kwargs):

  partition_name = kwargs["database_dir"]
  output_name = os.path.join( kwargs["working_dir"], os.path.basename(partition_name).split(".")[0]+".json")
  kwargs["output_name"] = output_name

  if kwargs.get("dask_worker_memory", "-1") == "-1":
    kwargs["dask_worker_memory"] = " "
    kwargs["queue_memory"] = ""
  else:
    kwargs["queue_memory"] = "--memory "+str( parse_memsize(kwargs["dask_worker_memory"])*kwargs["n_cpus"])
    kwargs["dask_worker_memory"] = " DASK_WORKER_MEMORY=%s "%kwargs["dask_worker_memory"]

  cmd_condor = 'python -m fragmenstein.external.condor_queue.send_to_condor --env_vars EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json ' \
        'PATH=/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin:$PATH ' \
        '%(dask_worker_memory)s N_CPUS=%(n_cpus)s DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT=30 --ncpus %(n_cpus)s %(queue_memory)s  "'

  cmd_args = ' -m fragmenstein.external.enamine_realDB_search.similarity_searcher_search_onePartition ' \
             '-d %(database_dir)s -o %(output_name)s  %(smilesFname)s'

  if "metric" in kwargs:
    cmd_args += " --metric %s "%kwargs["metric"]

  if "n_hits_per_smi" in kwargs:
    cmd_args += " --n_hits_per_smi %s"% kwargs["n_hits_per_smi"]

  if "backend" in kwargs:
    cmd_args += " --backend %s"% kwargs["backend"]

  if kwargs.get("verbose", False):
    cmd_args += " -v "

  if not run_locally:
    python = " /data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python "
    cmd = cmd_condor + python + cmd_args +'"'
  else:
    cmd = "python " + cmd_args


  cmd = cmd%kwargs

  print(cmd)
  check_call( cmd, shell=True)

def globalSearch():
  import sys
  import argparse
  from fragmenstein.utils.cmd_parser import ArgumentParser

  parser = ArgumentParser(prog="fast_similarity_search",
                          description="Find the K most similar compounds in the database")

  parser.add_argument('smiles_query', nargs=None, type=argparse.FileType('r'), default=sys.stdin,
                      help="smiles file with as many smiles as rows")

  parser.add_argument('-d', '--database_dir', help="the directory where compounds database was compiled", required=True)

  parser.add_argument('-m', '--metric', nargs=None, choices=["Tanimoto", "Traversky"], default="Traversky",
                      required=False,
                      help="metric to use")

  parser.add_argument('-n', '--n_hits_per_smi', type=int, default=30,
                      help="K highest score molecules to retrieve per query smiles ")

  parser.add_argument('-b', '--backend', choices=["numpy", "numba"], default="numba", required=False,
                      help="computation backend ")

  parser.add_argument('-w', '--working_dir', type=str, required=True,
                      help="The directory where per partition results will be saved")

  parser.add_argument('-l', '--run_locally', action="store_true",
                      help="run computations locally instead submitting to condor")

  parser.add_argument('-v', '--verbose', action="store_true", default=False,
                      help="Print to stdout working information ")

  args = parser.parse_args()
  query_smi_str = args.smiles_query.read()
  query_smi_list =query_smi_str.splitlines()
  assert len(query_smi_list) > 0, "Error, no smiles provided"

  kwargs = vars( args )

  if  not kwargs["run_locally"]:
    tmpdir = SHARED_TMP
  else:
    tmpdir = tempfile.tempdir

  with tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, suffix=".smi", delete=False) as f:
    f.write( query_smi_str )
    f.flush()
    kwargs["smilesFname"] = f.name

    del kwargs["smiles_query"]
    db_partitions = findDB_partitions(args.database_dir)

    for partitionDir in db_partitions:
      kwargs["database_dir"] = partitionDir
      launch_searcher( **kwargs )


if __name__ == "__main__":
  globalSearch()


'''
# Z2770753574
echo "CCNC(=O)C(F)(F)C(=O)NCC" | python -m fragmenstein.external.enamine_realDB_search.similarity_searcher_search_allPartitions -d ../../enamine/fingerprints  --n_cpus=4 -w ~/tmp/kkdir/ -

'''