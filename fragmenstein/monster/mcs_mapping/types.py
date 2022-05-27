from rdkit import Chem
from rdkit.Chem import rdFMCS
# TypedDict & Unpack fixed in .legacy:
from typing import Dict, List, Tuple, Optional, TypeVar, Sequence, TypedDict, Unpack  # noqa

IndexMap = TypeVar('IndexMap', bound=Sequence[Tuple[int, int]])
# this is not visible via help as IndexMap is technically an instance
IndexMap.__doc__ = """
Sequence (Tuple or List) of tuples of two indices, the first is the hit index, the second is the followup
"""


class BasicFMCSMode(TypedDict):
    """
    These are the acceptable types for the modes of FindMCS
    """
    atomCompare: rdFMCS.AtomCompare
    bondCompare: rdFMCS.BondCompare
    ringCompare: rdFMCS.RingCompare
    matchChiralTag: bool
    matchValences: bool
    completeRingsOnly: bool
    ringMatchesRingOnly: bool
    threshold: float
    maximizeBonds: bool
    timeout: int
    verbose: bool
    seedSmarts: str


class ExtendedFMCSMode(BasicFMCSMode):
    """
    These are the acceptable types for the modes of FindMCS, plus three that can be configured in rdFMCS.MCSParameters
    """
    matchStereo: bool
    matchFormalCharge: bool
    maxDistance: float
