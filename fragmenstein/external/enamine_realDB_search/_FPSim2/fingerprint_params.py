from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

FINGERPRINT_NBITS=2048
FINGERPRINT_RADIUS=2048
HD5_TEMPLATE="%s.fingerprints.h5"

N_LINES_PER_CHUNK = int(2e5) # int(1e6)
SPOOLING_TIME = 5