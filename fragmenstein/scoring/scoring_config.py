import os
from collections import OrderedDict
from datetime import date

N_WORKERS=4
N_THREADS_PER_WORKERS=1

MPRO_RAW_DATA_DIR = os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro")
RAW_SUBMISSIONS_TABLE_URL = "https://raw.githubusercontent.com/postera-ai/COVID_moonshot_submissions/master/covid_submissions_all_info.csv"
RAW_SUBMISSIONS_TABLE_PATH = os.path.expanduser("~/tmp/covid_submissions_all_info.csv")
MPRO_HITS_DIR = os.path.expanduser("~/oxford/tools/Fragmenstein/fragmenstein/mpro/data/hit_mols")

FEATURES_DIR = os.path.expanduser("~/oxford/myProjects/diamondCovid/scoring/computedFeatures")
FRAGMENSTEIN_OUT_DIR = os.path.expanduser("~/oxford/myProjects/diamondCovid/scoring/fragmenstein_out")

TMP_DIR = os.path.expanduser("~/tmp/")

SDF_SCORES_HEADER_INFO = OrderedDict([
    ("submitter_name", "Ruben Sanchez Garcia"),
    ("submitter_email", "ruben.sanchez-garcia@stats.ox.ac.uk"),
    ("submitter_institution", "University of Oxford"),
    ("generation_date", date.today().strftime("%Y-%m-%d")),
    ("method", "scoring_trial"),
    ("ref_url", "http://www.notprovided.com"),
    ("Name", "The name of the compound"),
    ("original SMILES", "The SMILES of the compound"),
    ("ref_pdb", "x0830_0"), # TODO: CURRENTLY, it works only with a ref_pdb for all molecules, the ref_pdb should be already available.
    ("score_xcos", "Score computed using Warren's XCOS") #TODO: autodetect
])