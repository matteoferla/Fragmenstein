import numba
import numpy as np


def jaccard_vectorized(x,y):
  intersections_count= x.astype(np.float32) @ y.T
  counts_x = np.sum(x, axis=-1)
  counts_y = np.sum(y, axis=-1)
  sums = counts_x.reshape(-1,1)+counts_y.reshape(-1,1).T
  union = sums - intersections_count
  return intersections_count / union


@numba.jit( nopython=True, cache=True)
def jaccard_numba(query_fp, db_fp):

    n11 = np.sum((query_fp == 1) & (db_fp == 1))
    n00 = np.sum((query_fp == 0) & (db_fp == 0))
    jac = n11 / (query_fp.shape[0] - n00)
    return jac


@numba.jit( nopython=True, cache=True)
def traversky_numba(query_fp, db_fp, alpha=0.3, beta=0.7):

    query_1 = (query_fp == 1)
    n11 = np.sum( query_1 & (db_fp == 1))
    n10 = np.sum( query_1 & (db_fp == 0))
    n01 = np.sum((query_fp == 0) & (db_fp == 1))
    res = n11 / (n11 + alpha*n10 +beta*n01)
    return res

def testTanimoto():
  from rdkit import DataStructs
  from fragmenstein.external.enamine_realDB_search.compute_fingerprints import get_fingerPrint

  smi1 = 'O=C(NCC(O)C1=CSC=C1)C1CCN(C2=CC=NC=C2)C1'
  smi2 = 'CCC(NC(=O)C1CCN(C2=CC=NC=C2)C1)C(=O)OC'

  fp1 = get_fingerPrint(smi1)
  fp2 = get_fingerPrint(smi2)

  sim = DataStructs.FingerprintSimilarity(fp1, fp2)
  print(sim)


def testTraversky():
  from rdkit import DataStructs
  from fragmenstein.external.enamine_realDB_search.compute_fingerprints import get_fingerPrint

  merge_smi = 'Cc1ccc([C@@](N)(O)OOCc2ccccc2)cc1F'
  f1_smi = 'CC=1C=CC(CS(=O)(=O)N)=CC1'
  f2_smi = 'OC=1C=CC(NC(=O)CCC=2C=CC=CC2)=CC1'

  fp_merge = get_fingerPrint(merge_smi)
  fp1 = get_fingerPrint(f1_smi)
  fp2 = get_fingerPrint(f2_smi)

  traversky_params = (0.3, 0.7)
  sim1 = DataStructs.FingerprintSimilarity(fp_merge, fp1,
                                           metric=lambda fp1, fp2: DataStructs.TverskySimilarity(fp1, fp2,
                                                                                                 *traversky_params))
  print(sim1)
  from fragmenstein.external.enamine_realDB_search.compute_fingerprints import get_fingerPrint_as_npBool
  print( traversky_numba(get_fingerPrint_as_npBool(merge_smi), get_fingerPrint_as_npBool(f1_smi), *traversky_params) )

  sim2 = DataStructs.FingerprintSimilarity(fp_merge, fp2,
                                           metric=lambda fp1, fp2: DataStructs.TverskySimilarity(fp1, fp2,
                                                                                                 *traversky_params))
  print(sim2)
  print( traversky_numba(get_fingerPrint_as_npBool(merge_smi), get_fingerPrint_as_npBool(f2_smi), *traversky_params) )

if __name__ == "__main__":
  testTraversky()