import os
import json
import dask.bag as db

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

  filenames = [os.path.join(path_to_jsons, fname) for fname in os.listdir(path_to_jsons)]
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


  parser.add_argument('-v', '--verbose', action="store_true", default=False,
                      help="Print to stdout working information ")

  args = parser.parse_args()

  final_search = combine_search_jsons(args.partitions_results_dir)
  print(final_search)
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