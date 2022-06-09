from ._no_blending import _MonsterNone
from typing import Dict, List, Tuple, Optional, Unpack, Set  # noqa: F401
from ..positional_mapping import GPM
import itertools
from ..mcs_mapping import IndexMap, ExtendedFMCSMode
from copy import deepcopy
from rdkit import Chem
from ..unmerge_mapper import Unmerge

class _MonsterExpand(_MonsterNone):
    """
    A variant of no_blend mode, with a focus on expansion
    """
    def by_expansion(self, primary_name:Optional[str]=None, min_mode_index:int=0) -> Chem.Mol:
        """
        Get the maps. Find the map with the most atoms covered.
        Use that map as the base map for the other maps.
        """
        # -------------- Get the primary hit ----------------------------
        primary_maps: List[Dict[int, int]]
        # primary_name as None chooses one, else the primary name provided is used:
        primary_name, primary_maps = self._get_primary_maps(primary_name)
        # -------------- Get the secondary hits --------------------------------
        # positional_overlap is used by ``_expand_hit_atom_map_by_overlap``
        # which is called by ``_get_unmerge_expansions``
        positional_overlaps: Dict[Tuple[str, str], Dict[int, int]] = self._compute_overlaps()
        unmergers: List[Unmerge] = self._get_unmerge_expansions(primary_name,
                                                                primary_maps,
                                                                positional_overlaps,
                                                                min_mode_index)
        # -------------- Sort the unmergers --------------------------------
        self.positioned_mol, self.mol_options = self._place_unmerger_expansions(unmergers)
        return self.positioned_mol

    def _get_primary_maps(self, primary_name:Optional[str]=None) -> Tuple[str, List[Dict[int, int]]]:
        """
        The primary hit is the hit will most in common with the placed molecule.

        :param primary_name:
        :return:
        """
        if primary_name is None:
            # the list is [{hit_atom_idx: template_atom_idx}, ...]
            maps: Dict[str, List[Dict[int, int]]] = self._compute_maps(broad=True)
            # get the largest maps (not the number of maps which would be `len(l)`)
            get_size = lambda l: len(l[0]) if len(l) else 0  # noqa: E731 Guido doesn't like lambda, but I do
            max_size = max(map(get_size, maps.values()))
            # sorted_maps: Dict[str, List[Dict[int, int]]] = dict(sorted(maps.items(),
            #                                                            key=lambda x: get_size(x[1]),
            #                                                            reverse=True))
            biggest_maps = {k: v for k, v in maps.items() if get_size(v) == max_size}
            # choose the first map
            primary_name = list(biggest_maps.items())[0][0]
            primary_maps: List[Dict[int, int]] = biggest_maps[primary_name]
        else:
            primary: Chem.Mol = self.get_hit_by_name(primary_name)
            primary_maps: List[Dict[int, int]] = self._compute_hit_maps(primary, broad=True)
        self.journal.debug(f"Primary hit: {primary_name} with {len(primary_maps)} Primary maps: {primary_maps}")
        return primary_name, primary_maps

    def _get_unmerge_expansions(self, primary_name, primary_maps,
                                positional_overlaps,
                                min_mode_index) -> List[Unmerge]:
        """
        Calls _perform_unmerge which calls Unmerge
        """
        unmergers = []
        # the no_blend mode does the unmerged based on a dict of optional maps,
        # i.e. the maps do not affect each other. Here it is important that they do.
        # hence each primary map is converted into a set of unmerge maps and the best wins.
        # if there is only one hit, then the primary map is the only unmerge map...
        if len(self.hits) == 1:
            return [self._perform_unmerge(maps={primary_name: primary_maps},
                                          n_poisonous=3,
                                          primary_name=primary_name
                                          )]
        # case: multiple hits
        for primary_map in primary_maps:  #: Dict[int, int]
            # iterate over the hit map and expand to all overlapping atoms
            self.journal.debug(f'primary_map: {primary_map}')
            exp_map: Dict[str, Dict[int, int]] = self._expand_hit_atom_map_by_overlap(primary_name,
                                                                                      primary_map,
                                                                                      positional_overlaps,
                                                                                      self.custom_map)
            self.journal.debug(f'initial expanded map (primary + overlaps): {exp_map}')
            exp_maps = {primary_name: [primary_map]}  # only one primary map!
            accounted_for: Set[int] = {i for i in primary_map.values() if i >= 0}
            # get the maps that are not the primary map
            for other in self.hits:
                other_name: str = other.GetProp('_Name')
                if other_name == primary_name:
                    continue
                mappings:List[Dict[int, int]]
                mode: ExtendedFMCSMode
                mappings, mode = self.get_mcs_mappings(other, self.initial_mol, min_mode_index, exp_map)
                # drop any that are redundant with the primary hit
                mappings = [d for d in mappings if len(set(d.values()) - accounted_for) > 0]
                exp_maps[other_name] = mappings
                self.journal.debug(f'candiate expanded maps: {exp_maps} following: {other_name}')
            # {h: f for h, f in .items() if h >= 0 and f >= 0}
            unmergers.append(self._perform_unmerge(maps=exp_maps,
                                                   n_poisonous=3,
                                                   primary_name=primary_name))
        return unmergers

    def _place_unmerger_expansions(self, unmergers: List[Unmerge]) -> Tuple[Chem.Mol, List[Chem.Mol]]:
        scores: List[int] = []
        mol_options = []
        best_mol = None
        for unmerger in unmergers:
            n_off_atoms: int = unmerger.offness(unmerger.combined, unmerger.combined_map)
            scores.append(len(unmerger.combined_map) - 3 * n_off_atoms)
        max_score = max(scores)
        # if they came out equal keep both...
        for score, unmerger in zip(scores, unmergers):
            if score != max_score:
                continue
            positioned_mol, inner_options = self._place_unmerger(unmerger)
            mol_options.extend(inner_options)
            if best_mol:  # it might get sorted again... so it is not important
                mol_options.insert(0, positioned_mol)
            else:
                best_mol = positioned_mol
        return best_mol, mol_options

    def _compute_overlaps(self) -> Dict[Tuple[str, str], Dict[int, int]]:
        positional_overlaps: Dict[Tuple[str, str], Dict[int, int]] = {}
        for mol1, mol2 in itertools.combinations(self.hits, 2):
            mol1_name:str = mol1.GetProp('_Name')
            mol2_name:str = mol2.GetProp('_Name')
            gpm = GPM.get_positional_mapping(mol1, mol2)
            positional_overlaps[(mol1_name, mol2_name)] = gpm
            positional_overlaps[(mol2_name, mol1_name)] = gpm
        return positional_overlaps


    def _include_missing_hits(self, custom_map: Dict[str, Dict[int, int]]) -> None:
        for hit in self.hits:
            name = hit.GetProp('_Name')
            if name not in custom_map:
                custom_map[name] = {}

    def _expand_hit_atom_map_by_overlap(self,
                                     hit_name: str,
                                     hit_atom_map: Dict[int, int],
                                     positional_overlaps: Dict[Tuple[str, str], Dict[int, int]],
                                     custom_map: Dict[str, Dict[int, int]]) -> Dict[str, Dict[int, int]]:
        """
        Expanded the custom_map by adding all atoms that are covered by the hit_atom_map.

        :param hit_name:
        :param hit_atom_map:
        :param positional_overlaps:
        :param custom_map:
        :param mode:
        :return: custom_map
        """
        expanded: Dict[str, Dict[int, int]] = deepcopy(custom_map)
        expanded[hit_name] = hit_atom_map
        self._include_missing_hits(expanded)
        for hit_atom_idx, template_atom_idx in hit_atom_map.items():
            for other in self.hits:
                other_name:str = other.GetProp('_Name')
                if other_name == hit_name:
                    continue
                # ------------- deal with atoms that overlap --------------------------
                overlaps: Dict[int, int] = positional_overlaps[(hit_name, other_name)]
                if hit_atom_idx in overlaps:
                    # ignore the special overrides
                    if template_atom_idx < 0:
                        continue
                    if hit_atom_idx < 0:
                        continue
                    if overlaps[hit_atom_idx] in expanded[other_name]:
                        continue
                    expanded[other_name][overlaps[hit_atom_idx]] = template_atom_idx
                # ------------- deal with atoms that do not overlap ------------------
                elif template_atom_idx in expanded[other_name].values():
                    pass  # there is a mapping already ?!
                else:  # damn the template_atom_idx
                    expanded[other_name][-2-hit_atom_idx] = template_atom_idx
        return expanded

