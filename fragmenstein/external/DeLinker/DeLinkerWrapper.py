import os, sys
import tempfile
import warnings

from rdkit_constrained_conformer_generation import gen_confs_constrained

from fragmenstein.external import ExternalToolImporter

warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore',category=FutureWarning)

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from matplotlib import pyplot as plt


from collections import defaultdict

# THIS_SCRIPT_PARENTDIR= os.path.dirname(__file__)
# sys.path.append( os.path.join(THIS_SCRIPT_PARENTDIR, "DeLinker_to_remove") )
# sys.path.append( os.path.join(THIS_SCRIPT_PARENTDIR, "DeLinker_to_remove", "analysis") )
# sys.path.append( os.path.join(THIS_SCRIPT_PARENTDIR, "DeLinker_to_remove", "examples") )
#
# DEEP_MODEL_FNAME= os.path.abspath(os.path.join(THIS_SCRIPT_PARENTDIR, "DeLinker_to_remove", "models/pretrained_DeLinker_model.pickle"))
#
# from DeLinker_test import DenseGGNNChemModel
# import frag_utils
# import rdkit_conf_parallel
# from data.prepare_data import read_file, preprocess
# import example_utils


(DeLinker_test, frag_utils, rdkit_conf_parallel,
                data_prepare_data, example_utils)=  ExternalToolImporter.import_tool("DeLinker",
                                                      ["DeLinker_test", "frag_utils", "rdkit_conf_parallel", "data.prepare_data", "example_utils"])

DEEP_MODEL_FNAME= os.path.abspath(os.path.join(ExternalToolImporter.get_rootdir("DeLinker"), "models/pretrained_DeLinker_model.pickle"))


def loadMol(fname, removeDummy=True, removeHs= True):
  mol = Chem.MolFromMolFile(fname)

  if removeDummy:
    edit_mol = Chem.RWMol(mol)
    for i in range(mol.GetNumAtoms()):
      if mol.GetAtomWithIdx(i).GetSymbol() == DUMMY_SYMBOL:
        edit_mol.RemoveAtom(i)

      mol= edit_mol.GetMol()

  if removeHs:
    mol= Chem.RemoveHs( mol )

  return mol

# How many cores for multiprocessing
n_cores = 6
# Whether to use GPU for generating molecules with DeLinker_to_remove
use_gpu = False
DUMMY_SYMBOL="*"

frag_1_path="/home/ruben/oxford/tools/Fragmenstein/fragmenstein/mpro/data/hit_mols/Mpro-x1458.mol"
frag_2_path="/home/ruben/oxford/tools/Fragmenstein/fragmenstein/mpro/data/hit_mols/Mpro-x0104.mol"

frag_1_sdf = loadMol(frag_1_path)
frag_1_smi = Chem.MolToSmiles(frag_1_sdf)

frag_2_sdf = loadMol(frag_2_path)
frag_2_smi = Chem.MolToSmiles(frag_2_sdf)

# print( Chem.MolToPDBBlock( frag_1_sdf ) )
# img = Draw.MolsToGridImage([Chem.MolFromSmiles(frag_1_smi), Chem.MolFromSmiles(frag_2_smi)], molsPerRow=2, subImgSize=(300, 300))
# plt.imshow(img); plt.show()

combo_no_exit = Chem.CombineMols(frag_1_sdf, frag_2_sdf)
combo = Chem.CombineMols(combo_no_exit, Chem.MolFromSmiles("*.*"))

# for i in range(combo.GetNumAtoms()):
#   atom= combo.GetAtomWithIdx(i)
#   print(i, atom.GetSymbol(), atom.GetExplicitValence() )

combo_2d = Chem.Mol(combo)
AllChem.Compute2DCoords(combo_2d)
mol_with_idxs= example_utils.mol_with_atom_index(combo_2d)
plt.imshow( Draw.MolsToGridImage([mol_with_idxs], molsPerRow=1) ); plt.show()


atom_idx_1 = 26
atom_idx_2 = 12


edcombo = Chem.EditableMol(combo)

num_heavy_atoms = combo.GetNumHeavyAtoms()

edcombo.AddBond(num_heavy_atoms, atom_idx_1, order=Chem.rdchem.BondType.SINGLE)
edcombo.AddBond(num_heavy_atoms+1, atom_idx_2, order=Chem.rdchem.BondType.SINGLE)
editedcombo = edcombo.GetMol()


AllChem.Compute2DCoords(editedcombo)



AllChem.Compute2DCoords(editedcombo)
Chem.SanitizeMol(editedcombo )

# plt.imshow( Draw.MolsToGridImage([editedcombo], molsPerRow=1) ); plt.show()


mol_to_link = edcombo.GetMol()
Chem.SanitizeMol(mol_to_link)

# Convert exit vectors to carbons for conformer generation
du = Chem.MolFromSmiles('*')
mol_to_link_carbon = AllChem.ReplaceSubstructs(mol_to_link,du,Chem.MolFromSmiles('C'),True)[0]
Chem.SanitizeMol(mol_to_link_carbon)
# Generate conformer
mol_to_link_carbon = Chem.AddHs(mol_to_link_carbon)
AllChem.ConstrainedEmbed(mol_to_link_carbon, combo_no_exit, randomseed=42)
mol_to_link_carbon = Chem.RemoveHs(mol_to_link_carbon)

# Add this conformer to the two unlinked fragments
conf = mol_to_link.GetConformer()
ref_conf = mol_to_link_carbon.GetConformer()
for i in range(mol_to_link_carbon.GetNumAtoms()):
    pos = list(ref_conf.GetAtomPosition(i))
    conf.SetAtomPosition(i, pos)
conf.SetId(0)
_ = mol_to_link.AddConformer(conf)

# plt.imshow( Draw.MolsToGridImage([mol_to_link], molsPerRow=1) ); plt.show()


# Get distance and angle between fragments
dist, ang = frag_utils.compute_distance_and_angle(mol_to_link, "", Chem.MolToSmiles(mol_to_link))

print(Chem.MolToSmiles(mol_to_link), dist, ang )

cwdir= os.getcwd()
with tempfile.TemporaryDirectory() as tmpdirname:
  os.chdir(tmpdirname)
  data_path = "./fragments_test_data.txt"
  with open(data_path, 'w') as f:
    f.write("%s %s %s" % (Chem.MolToSmiles(mol_to_link), dist, ang))
  raw_data = read_file(data_path)
  preprocess(raw_data, "zinc", "fragments_test", True)

  if not use_gpu:
      os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

  # Arguments for DeLinker_to_remove
  args = defaultdict(None)
  args['--dataset'] = 'zinc'
  args['--config'] = '{"generation": true, \
                       "batch_size": 1, \
                       "number_of_generation_per_valid": 50, \
                       "min_atoms": 2, "max_atoms": 6, \
                       "train_file": "molecules_fragments_test.json", \
                       "valid_file": "molecules_fragments_test.json", \
                       "output_name": "DeLinker_example_generation.smi"}'
  args['--freeze-graph-model'] = False
  args['--restore'] = DEEP_MODEL_FNAME

  # Setup model and generate molecules
  model = DeLinker_test.DenseGGNNChemModel(args)
  model.train()
  # Free up some memory
  del model

  # Load molecules
  generated_smiles = frag_utils.read_triples_file("./DeLinker_example_generation.smi")
os.chdir(cwdir)
in_mols = [smi[1] for smi in generated_smiles]
frag_mols = [smi[0] for smi in generated_smiles]
gen_mols = [smi[2] for smi in generated_smiles]

du = Chem.MolFromSmiles('*')
clean_frags = [Chem.MolToSmiles(Chem.RemoveHs(AllChem.ReplaceSubstructs(Chem.MolFromSmiles(smi),du,Chem.MolFromSmiles('[H]'),True)[0])) for smi in frag_mols]

print("Deep learning done!")


# Check valid
results = []
for in_mol, frag_mol, gen_mol, clean_frag in zip(in_mols, frag_mols, gen_mols, clean_frags):
    if len(Chem.MolFromSmiles(gen_mol).GetSubstructMatch(Chem.MolFromSmiles(clean_frag)))>0:
        results.append([in_mol, frag_mol, gen_mol, clean_frag])

print("Number of generated SMILES: \t%d" % len(generated_smiles))
print("Number of valid SMILES: \t%d" % len(results))
print("%% Valid: \t\t\t%.2f%%" % (len(results)/len(generated_smiles)*100))

import re

linkers = list(map(frag_utils.get_linker, [Chem.MolFromSmiles(m[2]) for m in results],
                   [Chem.MolFromSmiles(m[3]) for m in results], [m[1] for m in results]
                   ))
# Standardise linkers
for i, linker in enumerate(linkers):
  if linker == "":
    continue
  linker = Chem.MolFromSmiles(re.sub('[0-9]+\*', '*', linker))
  Chem.rdmolops.RemoveStereochemistry(linker)
  linkers[i] = Chem.MolStandardize.canonicalize_tautomer_smiles(Chem.MolToSmiles(linker))
# Update results
for i in range(len(results)):
  results[i].append(linkers[i])

print("Done")

# Create dictionary of results
results_dict = {}
for res in results:
    if res[0]+'.'+res[1] in results_dict: # Unique identifier - starting fragments and original molecule
        results_dict[res[0]+'.'+res[1]].append(tuple(res))
    else:
        results_dict[res[0]+'.'+res[1]] = [tuple(res)]
# Check uniqueness
print("Unique molecules: %.2f%%" % (frag_utils.unique(results_dict.values())*100))

# Check if molecules pass 2D filters
filters_2d = frag_utils.calc_filters_2d_dataset(results,
                                                pains_smarts_loc= os.path.abspath(os.path.join(THIS_SCRIPT_PARENTDIR,
                                                                                  "DeLinker_to_remove/analysis/wehi_pains.csv")),
                                                n_cores=n_cores)

results_filt = []
for res, filt in zip(results, filters_2d):
  if filt[0] and filt[1] and filt[2]:
    results_filt.append(res)

print("Pass all 2D filters: \t\t\t\t%.2f%%" % (len(results_filt) / len(results) * 100))
print("Valid and pass all 2D filters: \t\t\t%.2f%%" % (len(results_filt) / len(generated_smiles) * 100))
print(
  "Pass synthetic accessibility (SA) filter: \t%.2f%%" % (len([f for f in filters_2d if f[0]]) / len(filters_2d) * 100))
print("Pass ring aromaticity filter: \t\t\t%.2f%%" % (len([f for f in filters_2d if f[1]]) / len(filters_2d) * 100))
print(
  "Pass SA and ring filters: \t\t\t%.2f%%" % (len([f for f in filters_2d if f[0] and f[1]]) / len(filters_2d) * 100))
print("Pass PAINS filters: \t\t\t\t%.2f%%" % (len([f for f in filters_2d if f[2]]) / len(filters_2d) * 100))

# Get unique molecules
print("Number molecules passing 2D filters:\t\t%d" % len(results_filt))
results_filt_unique = example_utils.unique_mols(results_filt)
print("Number unique molecules passing 2D filters:\t%d" % len(results_filt_unique))


new_gen_mols= [ x[2] for x in results_filt_unique ]
args= dict(smiles_full_mols= [new_gen_mols] , sdffile= None, reference_confs= [mol_to_link], maxconfs=50, rms_threshold=0.35, energy=10,
          tdist=0.75, smi_frags= [frag_mols], numcores=n_cores, jpsettings=False, verbose=True )

from joblib import dump
dumpFname= os.path.join(THIS_SCRIPT_PARENTDIR, "args.pkl")
print( dumpFname )
dump(args, dumpFname)

confs = gen_confs_constrained ( [gen_mols] , None, [mol_to_link], maxconfs=10, rms_threshold=0.35, energy=10,
          tdist=0.75, smi_frags= [frag_mols], numcores=n_cores, jpsettings=False)


print(confs)

dumpFname= os.path.join(THIS_SCRIPT_PARENTDIR, "molecules_generated.pkl")
print( dumpFname )
dump(confs, dumpFname)