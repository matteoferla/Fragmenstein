import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs import cDataStructs

FINGERPRINT_NBITS=2048


def get_fingerPrint(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol == None:
        return None
    finprin= rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=FINGERPRINT_NBITS, useChirality=0, useBondTypes=1, useFeatures=0)
    return finprin

def computeFingerprintStr(smi):
    # print(smi)
    finPrint = get_fingerPrint(smi)
    if finPrint is None:
        return b""
    finStr = cDataStructs.BitVectToBinaryText(finPrint)
    return finStr


def get_fingerPrint_as_npBool(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol == None:
        return None
    finprin= rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=FINGERPRINT_NBITS, useChirality=0, useBondTypes=1, useFeatures=0)

    fp_num = np.zeros((0,), dtype=np.bool)
    DataStructs.ConvertToNumpyArray(finprin, fp_num)
    return fp_num


def computeFingerprint_np_Str(smi):
    # print(smi)
    finPrint = get_fingerPrint_as_npBool(smi)
    if finPrint is None:
        return b""
    # finStr = finPrint.tobytes()
    finStr = np.packbits(finPrint).tobytes()
    return finStr

def decompressFingerprint_npStr(fpr_np):
    return np.unpackbits(np.frombuffer(fpr_np, dtype=np.uint8)).astype(bool)
