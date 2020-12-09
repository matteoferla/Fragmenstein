import os, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
from rdkit.Chem import MolStandardize

import numpy as np

from itertools import product
from joblib import Parallel, delayed
import re
from collections import defaultdict

from IPython.display import clear_output
IPythonConsole.ipython_useSVG = True

THIS_SCRIPT_PARENTDIR= os.path.dirname(__file__)
sys.path.append( os.path.join(THIS_SCRIPT_PARENTDIR, "DeLinker") )
sys.path.append( os.path.join(THIS_SCRIPT_PARENTDIR, "DeLinker", "analysis") )
sys.path.append( os.path.join(THIS_SCRIPT_PARENTDIR, "DeLinker", "examples") )

from DeLinker_test import DenseGGNNChemModel
import frag_utils
import rdkit_conf_parallel
from data.prepare_data import read_file, preprocess
import example_utils



# How many cores for multiprocessing
n_cores = 4
# Whether to use GPU for generating molecules with DeLinker
use_gpu = False

frag_1_path="/home/ruben/oxford/tools/Fragmenstein/fragmenstein/mpro/data/hit_mols/Mpro-x1458.mol"
frag_2_path="/home/ruben/oxford/tools/Fragmenstein/fragmenstein/mpro/data/hit_mols/Mpro-x0104.mol"
frag_1_sdf = Chem.MolFromMolFile(frag_1_path)
frag_1_smi = Chem.MolToSmiles(frag_1_sdf)

frag_2_sdf = Chem.MolFromMolFile(frag_2_path)
frag_2_smi = Chem.MolToSmiles(frag_2_sdf)

print( Chem.MolToPDBBlock( frag_1_sdf ) )