# TypedDict & Unpack fixed in .legacy:
from typing import Dict, List, Tuple, Optional, TypeVar, Sequence, TypedDict, Unpack, overload  # noqa
from functools import singledispatch
from .types import IndexMap, BasicFMCSMode, ExtendedFMCSMode
from rdkit.Chem import rdFMCS



@overload
def flip_mapping(mapping: Sequence[Tuple[int, int]]) -> Sequence[Tuple[int, int]]:
    # basic IndexMap
    ...


@overload
def flip_mapping(mapping: Sequence[Sequence[Tuple[int, int]]]) -> Sequence[Sequence[Tuple[int, int]]]:
    # list of IndexMap options
    ...

@overload
def flip_mapping(mapping: Dict[int, int]) -> Dict[int, int]:
    ...

@overload
def flip_mapping(mapping: Dict[str, Dict[int, int]]) -> Dict[str, Dict[int, int]]:
    ...

@overload
def flip_mapping(mapping: Dict[str, Sequence[Tuple[int, int]]]) -> Dict[str, Sequence[Tuple[int, int]]]:
    ...

def flip_mapping(mapping):
    """
    This is a temporary fix for the fact that many route run on mappings that go hit to followup or viceversa
    Moving forward these will be replaced with a single hit to followup scheme, of type IndexMap
    """
    get_first_key = lambda d: list(d.keys())[0] if len(d) > 0 else None
    # empty
    if len(mapping) == 0:
        return mapping
    # Sequence[Tuple[int, int]] (basic)
    elif isinstance(mapping, Sequence) and isinstance(mapping[0], tuple) and isinstance(mapping[0][1], int):
        return [(b, a) for a, b in mapping]
    # Sequence[Sequence[Tuple[int, int]]]  (options
    elif isinstance(mapping, Sequence) and \
            isinstance(mapping[0], Sequence):
        return [flip_mapping(submapping) for submapping in mapping]
    elif isinstance(mapping, Sequence):  # boh
        raise NotImplementedError(f'{mapping}: list of ???')
    # Dict[int, int]
    elif isinstance(mapping, dict) and isinstance(get_first_key(mapping), int):
        return {b: a for a, b in mapping.items()}
    # Dict[str, Dict[int, int]]
    elif isinstance(mapping, dict) and isinstance(get_first_key(mapping), str):
        return {name: flip_mapping(submapping) for name, submapping in mapping.items()}
    # Dict[str, Sequence[Tuple[int, int]]
    elif isinstance(mapping, dict) and isinstance(get_first_key(mapping), str):
        return {name: flip_mapping(submapping) for name, submapping in mapping.items()}
    else:
        raise NotImplementedError(f'{mapping}: dictionary of {get_first_key(mapping)}')

# ----------------

def transmute_FindMCS_parameters(
        **mode: Unpack[ExtendedFMCSMode]) -> rdFMCS.MCSParameters:  # noqa lowercase not applicable
    """
    The function ``rdFMCS.FindMCS`` has two ways of being used.
    In one, a series of arguments are passed,
    in another a ``rdFMCS.MCSParameters`` object is passed (a wrapped C++ structure).
    Unfortunately, there does not seem to be a way to transmute the former into the other.

    Hence, this function

    The ``params.AtomTyper`` and ``params.BondTyper`` members
    can be either

    * the enum ``rdFMCS.AtomCompare`` or ``rdFMCS.BondCompare``
    * or a subclass of ``rdFMCS.MCSBondCompare`` or ``rdFMCS.MCSAtomCompare``

    The ``rdFMCS.RingCompare`` appears to be absorbed into the ``params.BondCompareParameters``
    member.
    """
    params = rdFMCS.MCSParameters()
    # three integers: https://github.com/rdkit/rdkit/blob/b208da471f8edc88e07c77ed7d7868649ac75100/Code/GraphMol/FMCS/FMCS.h
    # they are not rdFMCS.AtomCompare rdFMCS.BondCompare and rdFMCS.RingCompare enums?
    # atom parameters
    atomCompare: rdFMCS.AtomCompare = mode.get('atomCompare', rdFMCS.AtomCompare.CompareElements)
    params.AtomTyper = atomCompare  # int or a callable
    params.AtomCompareParameters.MatchIsotope = atomCompare == rdFMCS.AtomCompare.CompareIsotopes
    params.AtomCompareParameters.CompleteRingsOnly = mode.get('completeRingsOnly', False)
    params.AtomCompareParameters.MatchChiralTag = mode.get('matchChiralTag', False)
    params.AtomCompareParameters.MatchValences = mode.get('matchValences', False)
    params.AtomCompareParameters.RingMatchesRingOnly = mode.get('ringMatchesRingOnly', False)
    # bond parameters
    bondCompare: rdFMCS.BondCompare = mode.get('bondCompare', rdFMCS.BondCompare.CompareOrder)
    ringCompare: rdFMCS.RingCompare = mode.get('ringCompare', rdFMCS.RingCompare.IgnoreRingFusion)
    params.BondTyper = bondCompare
    params.BondCompareParameters.CompleteRingsOnly = mode.get('completeRingsOnly', False)
    params.BondCompareParameters.MatchFusedRings = ringCompare != rdFMCS.RingCompare.IgnoreRingFusion
    params.BondCompareParameters.MatchFusedRingsStrict = ringCompare == rdFMCS.RingCompare.StrictRingFusion
    params.BondCompareParameters.RingMatchesRingOnly = mode.get('ringMatchesRingOnly', False)
    params.Threshold = mode.get('threshold', 1.0)
    params.MaximizeBonds = mode.get('maximizeBonds', True)
    params.Timeout = mode.get('timeout', 3600)
    params.Verbose = mode.get('verbose', False)
    params.InitialSeed = mode.get('seedSmarts', '')
    # parameters with no equivalence (i.e. made up)
    params.BondCompareParameters.MatchStereo = mode.get('matchStereo', False)
    params.AtomCompareParameters.MatchFormalCharge = mode.get('matchFormalCharge', False)
    params.AtomCompareParameters.MaxDistance = mode.get('maxDistance', -1)
    # params.ProgressCallback Depracated
    return params
