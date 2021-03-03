
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit import DataStructs
import argparse
import sys

from rdkit.DataStructs import similarityFunctions

similarityFunctions = {  name:fun for name, fun, __ in similarityFunctions}

def search_smi_in_mol_list(mol_or_smiles, molId_mol_list, similarity_fun, n_jobs=1, verbose=False):
    '''
    Perforrms a similarity search against a set of molecules

    :param mol_or_smiles: query molecule
    :param molId_mol_list:
    :param n_jobs:
    :param verbose:
    :return:
    '''
    if not isinstance(mol_or_smiles, Chem.Mol):
        query_mol = Chem.MolFromSmiles(mol_or_smiles)
    else:
        query_mol = mol_or_smiles

    similarity_fun = similarityFunctions[similarity_fun]

    computeFingerprint = Chem.RDKFingerprint
    query_fingerprint = computeFingerprint(query_mol)
    scores = Parallel(n_jobs = n_jobs)( delayed( lambda mold_id, mol_i: (mold_id, DataStructs.FingerprintSimilarity(query_fingerprint,
                                                                                        computeFingerprint(mol_i), metric= similarity_fun) )
                                                    )(  mold_id, mol_i ) for mold_id, mol_i in molId_mol_list )
    scores = sorted(scores, key= lambda x: x[-1])
    if verbose:
        print( "\n".join( ( "%s\t%s"%(s_id, s) for s_id, s  in scores ) ) )
    return scores
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('query', nargs=None, help="the SMILES of the query protein")
    parser.add_argument('-d', '--database', nargs=None, type=argparse.FileType('r'), default=sys.stdin, required=True,
                        help="file text with as many rows are molecules. Each row has two fields: mol_id,smiles. The two fields"
                             "are comma separated")
    parser.add_argument('-m', '--metric', nargs=None, choices=list(similarityFunctions.keys()), default="Tanimoto", required=False,
                        help="metric to use")
    parser.add_argument('-t', '--n_jobs', nargs=None, type=int, default=1, help="Number of jobs for parallel execution")

    args = parser.parse_args()
    db_smiles = args.database.read().splitlines()

    db_smiles = (id_smiles.split(",") for id_smiles in db_smiles)
    database_mols = Parallel(n_jobs=args.n_jobs)( delayed( lambda  mol_id, smile: (mol_id, Chem.MolFromSmiles(smile)))( mol_id, smile ) for mol_id, smile in db_smiles )

    search_smi_in_mol_list(args.query, molId_mol_list= database_mols, similarity_fun= args.metric, n_jobs = args.n_jobs, verbose=True)

'''
This is an example


python -c 'import sys, os
dataDir=os.path.expanduser("~/oxford/myProjects/diamondCovid/data/nsp13/aligned")
for dirname in os.listdir(dataDir):
  frag_id = dirname.split("-")[1].split("_bound")[0]
  smiles_fname = os.path.join(dataDir, dirname, dirname+"_smiles.txt")
  if os.path.exists(smiles_fname):
    with open(smiles_fname) as f:
      smiles= f.read().strip()
    print("%s,%s"%(frag_id, smiles))' \
 | python fragmenstein/external/smilesSearcher/smiles_searcher.py  -d - "S=C1N(C(=NN1)CSC2=CC=CC=C2[N+]([O-])=O)CC=C"

'''