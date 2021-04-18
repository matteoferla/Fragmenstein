import os
import json
import dask.bag as db
import re

from collections import defaultdict

import itertools
from rdkit import Chem
from scipy.spatial.distance import pdist

from fragmenstein.utils.io_utils import apply_func_to_files


def combine_search_jsons(path_to_jsons):
  def combine_two_dict(d1, d2, ):
    new_d = {}
    assert len(d1) == len(d2), "Error, both dicts should be equaly-sized"
    for query_smi, matches_list1 in d1.items():
      matches_list2 = d2[query_smi]
      assert len(matches_list1) == len(matches_list2), "Error, both list should be equaly-sized"
      final_matches = []
      for match1, match2 in zip(matches_list1, matches_list2):
        similarity1 = match1[0]
        similarity2 = match2[0]
        if similarity1 > similarity2:
          final_matches.append(match1)
        else:
          final_matches.append(match2)
      new_d[query_smi] = final_matches
    return new_d

  def load_json(fname):
    with open(fname) as f:
      return json.load(f)

  filenames = [os.path.join(path_to_jsons, fname) for fname in os.listdir(path_to_jsons) if fname.endswith(".json") ]
  if len(filenames) >1:
    bag = db.from_sequence(filenames[1:]).map(load_json).fold(combine_two_dict, initial=load_json(filenames[0]))
    final_search = bag.compute()
  else:
    final_search = load_json(filenames[0])
  return final_search


def collectResults():

  from fragmenstein.utils.cmd_parser import ArgumentParser

  parser = ArgumentParser(prog="collect_similarity_results",
                          description="Collects the results obtained from similarity_searcher_search_allPartitions")

  parser.add_argument('partitions_results_dir', type=str,
                      help="The directory containing json files with results")

  parser.add_argument('-o', '--output_name', type=str, required=False,
                      help="The fname for a json file where search results will be stored")

  parser.add_argument('-f', '--fragments_dir', type=str, required=False,
                      help='The directory where smiles files  r".*?(x[\w-]+)_smiles\.txt" are contained for the used queries' )

  parser.add_argument('-v', '--verbose', action="store_true", default=False,
                      help="Print to stdout found results")

  parser.add_argument( '--for_fragmenstein_fname', type=str, required=False, default=None,
                      help="fname that can be readily used with fragmenstein protocols_placeFragmenstein ")


  args = parser.parse_args()

  final_search = combine_search_jsons(args.partitions_results_dir)

  # print(final_search )

  if args.fragments_dir:
    smiles_file_pattern =  r".*?(x[\w-]+)_smiles\.txt$"
    def readTxt(fname):
      with open(fname) as f:
        molId = re.match(smiles_file_pattern, os.path.split(fname)[-1]).group(1)
        return f.read().strip(), molId
    smiles_to_id = apply_func_to_files(args.fragments_dir, smiles_file_pattern, readTxt)

    mol_file_pattern =  r".*?(x[\w-]+)\.mol$"
    def computeCentroid(fname):
      print(os.path.split(fname)[-1])
      molId = re.match(mol_file_pattern, os.path.split(fname)[-1]).group(1)
      mol = Chem.MolFromMolFile(fname)
      return molId, mol.GetConformer().GetPositions().mean(axis=0)

    fragIds_centroids = apply_func_to_files(args.fragments_dir, mol_file_pattern, computeCentroid)
    fragIds_centroids = dict(fragIds_centroids)

    smiles_to_id_dict = defaultdict(list)
    for key, val in smiles_to_id:
      smiles_to_id_dict[key].append(val)
    smiles_to_id = dict( smiles_to_id_dict )

    new_dict = {}
    for query, resultsList in final_search.items():
      fragIds_combinations = itertools.product( * map(lambda x: smiles_to_id[x], query.split(",")))

      def checkThrDist(fragmentIds, distThr=15):
        centroids = [ fragIds_centroids[fragId] for fragId in fragmentIds ]
        distsMat = pdist(centroids)
        return max(distsMat) < distThr

      compatible_fragIds = filter( checkThrDist, fragIds_combinations)

      for queryAsFragIds in compatible_fragIds:
        queryAsFragIds = ",".join(queryAsFragIds)
        new_dict[queryAsFragIds] = resultsList
    final_search = new_dict

    if args.for_fragmenstein_fname:
      smis_set = set([])
      assert args.fragments_dir is not None, "Error, if for_fragmenstein_fname provided, fragment names are required"
      with open(args.for_fragmenstein_fname, "w") as f:
        for queryIds, dataList in final_search.items():
          for sim, cid, smi in dataList:
            mol = Chem.MolFromSmiles(smi)
            smi = Chem.MolToSmiles(mol)
            if smi not in smis_set:
              f.write("%s\t%s\n" % ( smi, queryIds))
              smis_set.add(smi)

  # if args.fragments_json:
  #   with open(args.fragments_json) as f:
  #     hit_to_smiles = json.load( f)
  #   new_dict = {}
  #   for query, resultsList in final_search.items():
  #     query = tuple( map(lambda x: hit_to_smiles[x], query.split(",") ))
  #     new_dict[query] = resultsList

  if args.output_name:
    with open(args.output_name, "w") as f:
      json.dump(final_search, f)

  if args.verbose:
    for query_smi, dataList in final_search.items():
      print("> Query: %s" % query_smi)
      for sim, cid, smi in dataList:
        print("%.4f\t%s\t%s" % (sim, cid, smi))


if __name__ == "__main__":
  collectResults()

'''

python -m  fragmenstein.external.enamine_realDB_search.similarity_searcher_collect_results -v ~/tmp/kkdir/

'''