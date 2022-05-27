"""
The ``mcs_mapping`` module contains the classes and functions for mapping with the MCS algorithm,
restricted to a mapping that must be obeyed.
"""


from .types import IndexMap, BasicFMCSMode, ExtendedFMCSMode
from .compare_atom import VanillaCompareAtoms, SpecialCompareAtoms
from .utils import flip_mapping, transmute_FindMCS_parameters