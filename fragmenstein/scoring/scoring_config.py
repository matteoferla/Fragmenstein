import os

N_WORKERS=4
N_THREADS_PER_WORKERS=1

MPRO_RAW_DATA_DIR = os.path.expanduser("~/oxford/myProjects/diamondCovid/data/Mpro")
RAW_SUBMISSIONS_TABLE_URL = "https://raw.githubusercontent.com/postera-ai/COVID_moonshot_submissions/master/covid_submissions_all_info.csv"
RAW_SUBMISSIONS_TABLE_PATH = os.path.expanduser("~/tmp/covid_submissions_all_info.csv")
MPRO_HITS_DIR = os.path.expanduser("~/oxford/tools/Fragmenstein/fragmenstein/mpro/data/hit_mols")

FEATURES_DIR = os.path.expanduser("~/oxford/myProjects/diamondCovid/scoring/computedFeatures")
FRAGMENSTEIN_OUT_DIR = os.path.expanduser("~/oxford/myProjects/diamondCovid/scoring/fragmenstein_out")

TMP_DIR = os.path.expanduser("~/tmp/")