from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

FINGERPRINT_NBITS=2048


def get_fingerPrint(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol == None:
        return None
    finprin= rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=FINGERPRINT_NBITS, useChirality=0, useBondTypes=1, useFeatures=0)
    return finprin


def open_compound_db(directory):
    raise NotImplementedError()
