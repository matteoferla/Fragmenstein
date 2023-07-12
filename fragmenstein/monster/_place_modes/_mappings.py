import warnings
from typing import Optional, Dict, List, Tuple, Set, Unpack, Union  # noqa cf. .legacy monkeypatch

from rdkit import Chem
from rdkit.Chem import rdFMCS

from .._merge import _MonsterMerge
from ..mcs_mapping import SpecialCompareAtoms, IndexMap, ExtendedFMCSMode, transmute_FindMCS_parameters


class _MonsterMap(_MonsterMerge):
    def get_mcs_mappings(self,
                         hit: Chem.Mol,
                         followup: Chem.Mol,
                         min_mode_index: int = 0,
                         custom_map: Optional[Dict[str, Dict[int, int]]] = None
                         ) -> Tuple[List[Dict[int, int]], ExtendedFMCSMode]:
        """
        This is a curious method. It does a strict MCS match.
        And then it uses laxer searches and finds the case where a lax search includes the strict search.

        :param hit: query molecule
        :param followup: target/ref molecule
        :param min_mode_index: the lowest index to try (opt. speed reasons)
        :param custom_map: is the user defined hit name to list of tuples of hit and followup index pairs
        :return: mappings and mode
        """
        if not custom_map:
            # Issue 42:
            # Dict[str, IndexMap] = Dict[str, Sequence[Tuple[int, int]]]
            # not Dict[str, Dict[int, int]]
            custom_map: Dict[str, Dict[int, int]] = self.custom_map
        custom_map = self.fix_custom_map(custom_map)  # Issue 42
        hit_name: str = hit.GetProp('_Name')
        # -------------- Most strict mapping -------------------------------
        # run from strictest to laxest to find the strictest mapping that encompasses the custom map
        inverted_strict_i = 0
        for strict_i, mode in enumerate(reversed(self.matching_modes + [self.strict_matching_mode]), start=-1):
            # required for limiting the next iterator
            inverted_strict_i = len(self.matching_modes) - strict_i
            # this calls `self._get_atom_maps` does the MCS search
            strict_maps: List[IndexMap] = self._get_atom_maps(hit=hit,
                                                              followup=followup,
                                                              custom_map=custom_map,
                                                              **mode)
            # there is the possibility that the strict mapping does not allow the
            # provided custom_map
            # if so only the provided custom_map hits are used.
            if hit_name in self.custom_map:
                # wanted hit indices
                wanted: List[int] = self._get_required_indices_for_map(self.custom_map[hit_name])
                strict_maps: List[Dict[int, int]] = self._validate_vs_custom(strict_maps, wanted)
                if len(strict_maps) != 0:
                    # these maps are valid
                    break  # from the reverse loop...
        else:
            self.journal.warning('Provided mapping is very unfavourable... using that along for expanding the search')
            # unexpected pycharm warning as list({1:1}.items()) does give [[(1,1)]]
            # single choice list in this case
            strict_maps: List[Dict[int, int]] = [custom_map.get(hit_name, {}), ]
        self.journal.debug(f'`get_mcs_mappings` strict_maps for {hit_name}: {strict_maps}')
        # -------------- Expand mapping -------------------------------
        # go from laxest to strictest until one matches the strict form...
        for i, mode in enumerate(self.matching_modes):
            if i < min_mode_index:
                continue
            if i == inverted_strict_i:
                break
            # heme is 70 or so atoms & Any-Any matching gets stuck. So guestimate to avoid that:
            if hit.GetNumAtoms() > (50 + i * 10) or followup.GetNumAtoms() > (50 + i * 10):
                continue

            lax: List[IndexMap] = []
            for strict_map in strict_maps:  #: Dict[int, int]
                # `_get_atom_maps` does the MCS search constrained by `expanded_custom_map`
                expanded_custom_map: Dict[str, Dict[int, int]] = self.expand_custom_map(custom_map,
                                                                                        {hit_name: strict_map})
                lax.extend(self._get_atom_maps(hit=hit,
                                               followup=followup,
                                               custom_map=expanded_custom_map,
                                               **mode)
                           )
            if len(lax) == 0:
                continue
            else:
                return [dict(n) for n in lax], mode
        # The strict will have to do.
        return [dict(n) for n in strict_maps], self.strict_matching_mode  # tuple to dict

    def get_mcs_mapping(self, hit, followup, min_mode_index: int = 0) -> Tuple[Dict[int, int], dict]:
        """
        This is a weird method. It does a strict MCS match.
        And then it uses laxer searches and finds the case where a lax search includes the strict search.

        :param hit: query molecule
        :param followup: target/ref molecule
        :param min_mode_index: the lowest index to try (opt. speed reasons)
        :return: mapping and mode
        """
        ms, mode = self.get_mcs_mappings(hit, followup, min_mode_index)
        return ms[0], mode

    def _get_atom_maps(self,
                       hit: Chem.Mol,
                       followup: Chem.Mol,
                       custom_map: Optional[Dict[str, IndexMap]]=None,
                       **mode: Unpack[ExtendedFMCSMode]) -> List[Dict[int, int]]:

        """
        The ``mode`` are FindMCS arguments, but this transmutes them into parameter scheme

        The old method is now ``_get_atom_maps_OLD``.
        """
        if custom_map is None:
            custom_map: Dict[str, Dict[int, int]] = self.custom_map
        custom_map = self.fix_custom_map(custom_map)  # Issue 42
        parameters: rdFMCS.MCSParameters = transmute_FindMCS_parameters(**mode)
        # this looks odd, because the default parameters.AtomTyper is a atomcompare enum
        # and can be overridden by a callable class instance (of MCSAtomCompare)
        parameters.AtomTyper = SpecialCompareAtoms(comparison=parameters.AtomTyper,
                                                   custom_map=custom_map)
        res: rdFMCS.MCSResult = rdFMCS.FindMCS([hit, followup], parameters)
        matches: List[Dict[int, int]] = parameters.AtomTyper.get_valid_matches(parameters.AtomCompareParameters,
                                                                         common=Chem.MolFromSmarts(res.smartsString),
                                                                         hit=hit,
                                                                         followup=followup
                                                                         )
        return [dict(m) for m in matches]

    def _get_atom_maps_OLD(self, molA, molB, **mode: Unpack[ExtendedFMCSMode]) -> Set[Tuple[Tuple[int, int]]]:
        """
        The ``mode`` are FindMCS arguments.
        """
        warnings.warn('`_get_atom_maps_OLD` is no longer used.', category=DeprecationWarning)
        mcs = rdFMCS.FindMCS([molA, molB], **mode)
        common = Chem.MolFromSmarts(mcs.smartsString)
        matches = []
        # prevent a dummy to match a non-dummy, which can happen when the mode is super lax.
        is_dummy = lambda mol, at: mol.GetAtomWithIdx(at).GetSymbol() == '*'  # noqa: E731 is stupid

        def all_bar_dummy(Aat, Bat) -> bool:
            """it is okay to match a dummy with dummy only"""
            a_dummy: bool = is_dummy(molA, Aat)
            b_dummy: bool = is_dummy(molB, Bat)
            # xor
            return (a_dummy and b_dummy) or not (a_dummy or b_dummy)

        # prevent matching hydrogens
        is_hydrogen = lambda mol, at: mol.GetAtomWithIdx(at).GetSymbol() == 'H'  # noqa: E731 is stupid
        for molA_match in molA.GetSubstructMatches(common, uniquify=False):
            for molB_match in molB.GetSubstructMatches(common, uniquify=False):
                matches.append([(int(molA_at), int(molB_at)) for molA_at, molB_at in zip(molA_match, molB_match) if
                                all_bar_dummy(molA_at, molB_at) and
                                not is_hydrogen(molA, molA_at) and
                                not is_hydrogen(molB, molB_at)])
        # you can map two toluenes 4 ways, but two are repeats.
        matches = set([tuple(sorted(m, key=lambda i: i[0])) for m in matches])
        return matches

    def _get_atom_map(self, molA, molB, **mode) -> List[Tuple[int, int]]:
        return self._get_atom_maps(molA, molB, **mode)[0]

    def expand_custom_map(self,
                          custom_map: Dict[str, Dict[int, int]],
                          addend: Dict[str, Dict[int, int]]) \
            -> Dict[str, Dict[int, int]]:
        # not a IndexMap
        new_map = custom_map.copy()
        for k, add_mapping in addend.items():
            if k not in new_map:
                new_map[k] = dict(add_mapping)
            else:
                # custom_map gets priority
                new_map[k] = {**dict(add_mapping), **dict(custom_map[k])}
        return new_map

    def _validate_vs_custom(self, maps: List[Dict[int, int]], wanted_idx: List[int]) -> List[Dict[int, int]]:
        """return only the IndexMaps with all wanted_idx (wanted in the hit)
        This is not part of ``SpecialCompareAtoms.get_valid_matches``
        because there is a difference between must match and must be present as discussed in that method.
        """
        get_set = lambda mapping: set(mapping.keys()) if isinstance(mapping, dict) else {i for i, j in mapping}
        return [mapping for mapping in maps if len(set(wanted_idx) - get_set(mapping)) == 0]

    def _get_required_indices_for_map(self, custom_hit_map: Dict[int, int]) -> List[int]:
        """
        Returns the hit indices that must be present in the map.
        """
        return [h for h, f in custom_hit_map.items() if h >= 0 and f >= 0]
