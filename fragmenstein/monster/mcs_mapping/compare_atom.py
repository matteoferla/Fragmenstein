from rdkit import Chem
from rdkit.Chem import rdFMCS
# TypedDict & Unpack fixed in .legacy:
from typing import Dict, List, Union, Tuple, Optional, TypeVar, Sequence, TypedDict, Unpack  # noqa
from functools import singledispatchmethod  # monkeypatched by .legacy (for Py3.7)
import itertools

from .types import IndexMap, BasicFMCSMode, ExtendedFMCSMode


class VanillaCompareAtoms(rdFMCS.MCSAtomCompare):
    """
    Give than the following does not work, one has to do it in full:

    .. code-block:: python
        :caption: This will raise an error:
        super().__call__(parameters, mol1, atom_idx1, mol2, atom_idx2)

    This class replicates the vanilla functionality
    """

    def __init__(self, comparison: rdFMCS.AtomCompare = rdFMCS.AtomCompare.CompareAnyHeavyAtom):
        """
        Whereas the atomCompare is an enum, this is a callable class.
        But in parameters there is no compareElement booleans etc. only Isotope...
        In https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/FMCS/Wrap/testFMCS.py
        it is clear one needs to make one's own.
        """
        super().__init__()  # noqa the p_obect (python object aka ctypes.py_object) is not used
        self.comparison = comparison

    def __call__(self,  # noqa signature matches... it is just Boost being Boost
                 parameters: rdFMCS.MCSAtomCompareParameters,
                 mol1: Chem.Mol,
                 atom_idx1: int,
                 mol2: Chem.Mol,
                 atom_idx2: int) -> bool:
        a1: Chem.Atom = mol1.GetAtomWithIdx(atom_idx1)
        a2: Chem.Atom = mol2.GetAtomWithIdx(atom_idx2)
        # ------- isotope ------------------------
        if parameters.MatchIsotope and a1.GetIsotope() != a2.GetIsotope():  # noqa
            return False
        if sum([a1.GetAtomicNum() == 0, a2.GetAtomicNum() == 0]) == 1:
            return False
        elif self.comparison == rdFMCS.AtomCompare.CompareIsotopes and a1.GetIsotope() != a2.GetIsotope():  # noqa
            return False
        elif self.comparison == rdFMCS.AtomCompare.CompareElements and a1.GetAtomicNum() != a2.GetAtomicNum():  # noqa
            return False
        elif self.comparison == rdFMCS.AtomCompare.CompareAnyHeavyAtom \
                and (a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1):  # noqa
            return False
        elif self.comparison == rdFMCS.AtomCompare.CompareAny:
            pass
        # ------- valence ------------------------
        if parameters.MatchValences and a1.GetTotalValence() != a2.GetTotalValence():  # noqa
            return False
        # ------- chiral ------------------------
        if parameters.MatchChiralTag and not self.CheckAtomChirality(parameters, mol1, atom_idx1, mol2, atom_idx2):
            return False
        # ------- formal ------------------------
        if parameters.MatchFormalCharge and not self.CheckAtomCharge(parameters, mol1, atom_idx1, mol2, atom_idx2):
            return False
        # ------- ring ------------------------
        if parameters.RingMatchesRingOnly and not self.CheckAtomRingMatch(parameters, mol1, atom_idx1, mol2, atom_idx2):
            return False
        # ------- complete ------------------------
        if parameters.CompleteRingsOnly:
            # don't know its intended to be used
            pass
        # ------- distance ------------------------
        if parameters.MaxDistance:  # integer...
            # todo fill!
            pass
        return True


class SpecialCompareAtoms(VanillaCompareAtoms):
    """
    This works like the ``_get_atom_maps`` did prior to Fragmentein version 0.9.
    The mapping as discussed in GitHub issue #23 is in the format

        mapping = { 'hit1': {1:1,2:5} 'hit2': {3:3,4:4,4:6}}

    The hit index is first, followup index is the second.
    The index `-1` for a followup index is the same as not providing the hit index,
    it is described here solely for clarity not for use.

        mapping = { 'hit1': {1:1,2:5, 3:-1} 'hit2': {3:3,4:4,4:6}}

    The index `-2` for a followup index will result in the hit atom index
    not matching any followup index.

        mapping = { 'hit1': {1:1,2:5, 3:-2} 'hit2': {3:3,4:4,4:6}}

    If ``exclusive_mapping`` argument of __init__ is True,
    then if a followup index is present in one hit, but not in a second hit,
    then no atom of the second hit will match that followup atom.
    A negative index for a hit atom index means that no atom in that hit will match the
    corresponding followup index.

        mapping = { 'hit1': {1:1,2:5,-1:3, -2: 7} 'hit2': {3:3,4:4,4:6}}

    However, a positive integer on a different hit overrides it, therefore,
    in the above followup atom 3 cannot be matched to any atom in hit1, but will match
    atom 3 in hit2. Followup atom 7 will not match either.

    .. code-block:: python

        SpecialCompareAtoms(custom_map=mapping, exclusive_mapping=True)
    """

    def __init__(self,
                 comparison: rdFMCS.AtomCompare = rdFMCS.AtomCompare.CompareAnyHeavyAtom,
                 custom_map: Optional[Dict[str, Dict[int, int]]] = None, exclusive_mapping: bool = True):
        super().__init__(comparison=comparison)  # what is p_object?
        self.custom_map = self.fix_custom_map(custom_map)  # custom as in user-defined
        self.banned = self._get_strict_banned() if exclusive_mapping else self._get_lax_banned()

    def _get_strict_banned(self) -> List[int]:
        """
        A list of followup indices that cannot be unmapped
        called if exclusive_mapping is True
        """
        return [foll_idx for mapping in self.custom_map.values() for foll_idx in mapping.values()]

    def _get_lax_flipped_map(self) -> List[int]:
        """
        A list of followup indices that cannot be mapped as per negative hit index
        called if exclusive_mapping is False
        """
        return [foll_idx for mapping in self.custom_map.values() for hit_idx, foll_idx in mapping.items()
                if hit_idx < 0]

    def get_custom(self, hit_mol: Chem.Mol, hit_atom_idx: int) -> int:
        """
        What idx of followup corresponds to the ``hit_atom_idx`` index of ``hit_mol``?
        If nothing, -1 is returned
        """
        name: str = hit_mol.GetProp('_Name')
        return self.custom_map[name].get(hit_atom_idx, -1)

    def __call__(self,
                 parameters: rdFMCS.MCSAtomCompareParameters,
                 hit: Chem.Mol,
                 hit_atom_idx: int,
                 followup: Chem.Mol,
                 followup_atom_idx: int) -> bool:
        if followup.GetProp('_Name') == hit.GetProp('_Name'):
            # self to self mapping
            return super().__call__(parameters, hit, hit_atom_idx, followup, followup_atom_idx)
        if hit.GetProp('_Name') not in self.custom_map:
            # self to self mapping
            return super().__call__(parameters, hit, hit_atom_idx, followup, followup_atom_idx)
        hit_atom = hit.GetAtomWithIdx(hit_atom_idx)
        followup_atom = followup.GetAtomWithIdx(followup_atom_idx)
        symbols = {hit_atom.GetSymbol(), followup_atom.GetSymbol()}
        # ------- Custom -----------------------
        # get the user defined map target -- see docstring for bla bla
        custom: int = self.get_custom(hit, hit_atom_idx)
        if custom == followup_atom_idx:
            # it is the custom map!
            return True
        elif custom != -1:
            # a different index was given
            # or a non -1 negative number
            return False
        elif followup_atom_idx in self.banned:
            # banned is fully filled if exclusive_mapping was true during init
            # otherwise its user provided negatives
            return False
        else:
            # followup index not assigned
            pass
        # ------- Dummy ------------------------
        # dummy atom cannot match non-dummy atom:
        if '*' in symbols and len(symbols) > 1:
            return False
        # ------- protons ------------------------
        # proton cannot match non-proton:
        if 'H' in symbols and len(symbols) > 1:
            return False
        # ------- vanilla ------------------------
        return super().__call__(parameters, hit, hit_atom_idx, followup, followup_atom_idx)

    @singledispatchmethod
    def get_valid_matches(self,
                          parameters: rdFMCS.MCSAtomCompareParameters,
                          common: Chem.Mol,
                          hit: Chem.Mol,
                          followup: Chem.Mol) -> List[Dict[int, int]]:
        """
        Returns a list of possible matches, each a dictionary of hit to follow indices,
        that obey the criteria of the atomic comparison.
        (Formerly this was a IndexAtom)

        This however does not check if all the atoms in custom are present.
        For that, ``Monster._validate_vs_custom`` is called in the method ``Monster.get_mcs_mappings``,
        after calling ``Monster._get_atom_maps``, which calls this method.
        (``Monster.get_mcs_mappings`` tries different matching schema,
        while ``Monster._get_atom_maps`` is for one single scheme).
        The primary reason why this is so, is that there are two tiers of requirements:

        1. The custom map must be present vs
        2. The custom map may be present, but has to be obeyed.

        parameters can be rdFMCS.MCSAtomCompareParameters or rdFMCS.MCSParameters
        """

        matches = []
        # find the common matches and make sure they match each other:
        # GetSubstructMatches is independent of AtomCompare.
        for hit_match, followup_match in itertools.product(hit.GetSubstructMatches(common, uniquify=False),
                                                           followup.GetSubstructMatches(common, uniquify=False)):
            # re `map(int, hit_match)` I do not know under what condition is it not an int...
            # but it was so in previous iteration
            for h, f in zip(map(int, hit_match), map(int, followup_match)):
                if not self(parameters=parameters,
                            hit=hit,
                            followup=followup,
                            hit_atom_idx=h,
                            followup_atom_idx=f):
                    break  # one index does not match. None should match.
            else:
                matches.append(dict(zip(map(int, hit_match), map(int, followup_match))))
        # remove duplicates
        hash_dict = lambda match: hash(tuple(sorted(match.items(), key=lambda i: i[0])))
        unique_matches = {hash_dict(match): match for match in matches}
        return list(unique_matches.values())

    @get_valid_matches.register
    def _(self,
          parameters: rdFMCS.MCSParameters,
          common: Chem.Mol,
          hit: Chem.Mol,
          followup: Chem.Mol) -> List[IndexMap]:
        return self.get_valid_matches(parameters.AtomCompareParameters, common, hit, followup)

    def fix_custom_map(self,
                          custom_map: Dict[str, Union[Sequence[Tuple[int, int]], Dict[int, int]]]) \
            -> Dict[str, Dict[int, int]]:
        """
        This will be deprecated in the future. As Monster.fix_custom_map is better.

        Make sure its Dict[str, Dict[int, int]]

        There is a bit of confusion about the custom map.
        Converts the custom map from dict of lists of 2-element tuples to dict of dicts.
        """
        if custom_map is None:
            # in Monster {h.GetProp('_Name'): {} for h in self.hits}
            return {}
        assert isinstance(custom_map, dict), 'User defined map has to be mol name to Dict[int, int]'
        for name, hit_map in custom_map.items():
            custom_map[name] = dict(hit_map)
        return custom_map
    
