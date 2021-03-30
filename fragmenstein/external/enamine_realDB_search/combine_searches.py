import json
import os

import numpy as np
from dask import bag as db


def combine_two_chunk_searches(cs1, cs2):
    '''

    :param cs1: tuple( similarities_1, ids_1)
    :param cs2: tuple( similarities_2, ids_2)
    :return:
    '''
    assert cs1[0].shape == cs2[0].shape, "Error, union for different shapes not implemented"

    sim_concat =  np.concatenate([ cs1[0], cs2[0]], axis=1 )
    idxs_concat = np.concatenate([ cs1[1], cs2[1]], axis=1 )

    to_pick_idxs = np.argsort(-sim_concat, axis = -1)[:, :cs1[0].shape[1]]
    new_sim =  -1 * np.ones_like(cs1[0])
    new_idxs = -1 * np.ones_like(cs1[1])
    for i in range(sim_concat.shape[0]):
        new_sim[i,:] = sim_concat[i, to_pick_idxs[i, :]]
        new_idxs[i,...] = idxs_concat[i, to_pick_idxs[i, :], :]

    # print("s1"); print(cs1[0])
    # print("s2"); print(cs2[0])
    # print("res="); print(new_sim)

    return new_sim, new_idxs


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