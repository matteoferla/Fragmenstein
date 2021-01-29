
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit import DataStructs
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('query', nargs=None, help="the SMILES of the query protein")
parser.add_argument('-d', '--database', nargs=None, type=argparse.FileType('r'), default=sys.stdin, required = True)
parser.add_argument('-t', '--n_jobs', nargs=None, type=int, default=1, help="Number of jobs for parallel execution")

args = parser.parse_args()

db_smiles = args.database.read().splitlines()

db_smiles = ( id_smiles.split(",") for id_smiles in db_smiles )

query_mol = Chem.MolFromSmiles(args.query)

database_mols = Parallel(n_jobs=args.n_jobs)( delayed( lambda  mol_id, smile: (mol_id, Chem.MolFromSmiles(smile)))( mol_id, smile ) for mol_id, smile in db_smiles )


computeFingerprint = Chem.RDKFingerprint
query_fingerprint = computeFingerprint(query_mol)
scores = Parallel(n_jobs=args.n_jobs)( delayed( lambda mold_id, mol_i: (mold_id, DataStructs.FingerprintSimilarity(query_fingerprint,
                                                                                                computeFingerprint(mol_i)) )
                                                )(  mold_id, mol_i ) for mold_id, mol_i in database_mols )
scores = sorted(scores, key= lambda x: x[-1])
print( "\n".join( ( "%s\t%s"%(s_id, s) for s_id, s  in scores ) ) )


'''

python -c 'import sys, os
dataDir=os.path.expanduser("~/oxford/myProjects/diamondCovid/data/nsp13/aligned")
for dirname in os.listdir(dataDir):
  frag_id = dirname.split("-")[1].split("_bound")[0]
  smiles_fname = os.path.join(dataDir, dirname, dirname+"_smiles.txt")
  if os.path.exists(smiles_fname):
    with open(smiles_fname) as f:
      smiles= f.read().strip()
    print("%s,%s"%(frag_id, smiles))' \
 | python /home/ruben/oxford/tools/Fragmenstein/fragmenstein/external/smilesSearcher/smiles_searcher.py  -d - "S=C1N(C(=NN1)CSC2=CC=CC=C2[N+]([O-])=O)CC=C"

'''