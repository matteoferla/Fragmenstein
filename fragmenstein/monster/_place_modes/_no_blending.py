import operator
from typing import Optional, Dict, List, Tuple, Set, Sequence, Unpack, Union, Any, Iterable  # noqa cf. .legacy monkeypatch

from rdkit import Chem
from rdkit.Chem import AllChem

from ._refine import _MonsterRefine
from ..mcs_mapping import IndexMap, flip_mapping
from ..unmerge_mapper import Unmerge
from collections import Counter
from ...error import DistanceError, PoisonError

class _MonsterNone(_MonsterRefine):

    def no_blending(self, broad=False) -> None:
        """
        no merging is done. The hits are mapped individually. Not great for small fragments.
        """
        maps: Dict[str, List[Dict[int, int]]] = self._compute_maps(broad)
        # there is no primary hit in no_blending mode.
        unmerger = self._perform_unmerge(maps, n_poisonous=3 if broad else 0)
        self.unmatched = [m.GetProp('_Name') for m in unmerger.disregarded]
        self.positioned_mol, self.mol_options = self._place_unmerger(unmerger)


    def _perform_unmerge(self,
                         maps: Dict[str, List[Dict[int, int]]],
                         n_poisonous:int,
                         primary_name:Optional[str]=None) -> Unmerge:
        """
        The second third of the no_blending method.
        But also used by expansion mapping.

        :param maps:
        :return:
        """
        positive: Dict[str, List[Dict[int, int]]] = self._remove_negatives(maps)
        flipped_maps: Dict[str, List[Dict[int, int]]] =\
            {name: [flip_mapping(hm) for hm in hit_mappings] for name, hit_mappings in maps.items()}
        # todo flip round Unmerge. Unmerge wants name to dict of followup to hit... which is backwards.
        unmerger = Unmerge(followup=self.initial_mol,
                             mols=self.hits,
                             maps=flipped_maps,
                             no_discard=False)  # self.throw_on_discard?
        unmatched = [m.GetProp('_Name') for m in unmerger.disregarded]
        if n_poisonous and unmatched:
            self._remove_poisonous(flipped_maps, unmerger.poisonous_indices, n_poisonous, primary_name)
            retried_unmerger = Unmerge(followup=self.initial_mol,
                               mols=self.hits,
                               maps=flipped_maps,
                               no_discard=False)  # self.throw_on_discard?
            retried_unmatched = [m.GetProp('_Name') for m in retried_unmerger.disregarded]
            if len(retried_unmatched) < len(unmatched):
                self.journal.debug(f'Retried unmerge yield the matching of ' +
                                   f'{set(unmatched) - set(retried_unmatched)} ' +
                                   f'via the removal of {unmerger.poisonous_indices}')
                unmerger = retried_unmerger
                unmatched = retried_unmatched
        if not self.throw_on_discard:
            pass
        elif len(unmatched) == 0:
            pass
        elif len(unmerger.poisonous_indices):
            raise PoisonError(mol=unmatched, indices=unmerger.poisonous_indices)
        else:
            raise DistanceError(hits=unmatched)
        self.journal.debug(f'followup to scaffold {unmerger.combined_map}')
        return unmerger

    def _remove_poisonous(self,
                          flipped_maps,
                          poisonous_indices: List[int],
                          n_poisonous: int,
                          primary_name:Optional[str]=None) -> None:
        """
        By poisonous it is intended the few atoms that cause a mapping for placement to go wrong.

        :param flipped_maps: modified in place
        :param poisonous_indices: a list with repeated indices as unmerge.poisonous_indices: for commonness
        :param n_poisonous: max indices to remove
        :param primary_name: remove poisonous atoms from primary hit only
        :return:
        """
        for i, c in Counter(poisonous_indices).most_common(n_poisonous):  # followup index
            for name, followup2hits in flipped_maps.items():
                if primary_name is not None and name != primary_name:
                    continue
                for followup2hit in followup2hits:
                    if i not in followup2hit.keys():
                        continue
                    if name in self.custom_map and i in self.custom_map[name].values():
                        continue
                    self.journal.debug(f'Removing poisonous mapping: ' +
                                       f'hit {name} index {followup2hit[i]} to followup index {i}')
                    del followup2hit[i]

    def _place_unmerger(self, unmerger: Unmerge) -> Tuple[Chem.Mol, List[Chem.Mol]]:
        # ------------------ places the atoms with known mapping ------------------
        placed = self.place_from_map(target_mol=self.initial_mol,
                                     template_mol=unmerger.combined_bonded,
                                     atom_map=unmerger.combined_map,
                                     random_seed=self.random_seed)

        self.keep_copy(unmerger.combined, 'scaffold')
        self.keep_copy(unmerger.combined_bonded, 'chimera')

        alts = zip(unmerger.combined_bonded_alternatives, unmerger.combined_map_alternatives)
        placed_options = [self.place_from_map(target_mol=self.initial_mol,
                                             template_mol=mol,
                                             atom_map=mappa,
                                             random_seed=self.random_seed) for mol, mappa in alts]
        # ------------------ .posthoc_refine averages the overlapping atoms ------------------
        if unmerger.pick != -1:  # override present!
            self.journal.debug(f'override present: {unmerger.pick}')
            positioned_mol = self.posthoc_refine(placed)
            mol_options = [self.posthoc_refine(mol) for mol in placed_options]
        else:
            mols: List[Chem.Mol] = [self.posthoc_refine(mol) for mol in [placed] + placed_options]
            positioned_mol: Chem.Mol = self.get_best_scoring(mols)
            mols.remove(positioned_mol)
            mol_options: List[Chem.Mol] = mols
        return positioned_mol, mol_options

    def _compute_maps(self, broad: bool) -> Dict[str, List[Dict[int, int]]]:
        """
        Compute the mapping for each hit and returns a dict of name to list of mappings
        in the scheme of ``[{hit_atom_idx: template_atom_idx}, ...]``

        :param broad:
        :return:
        """
        return {hit.GetProp('_Name'): self._compute_hit_maps(hit, broad) for hit in self.hits}

    def _compute_hit_maps(self, template: Chem.Mol, broad: bool) -> List[Dict[int, int]]:
        """
        Calcualte the list of options of possible maps of hit
        This is called in the expansion tactic in ``_get_primary_maps``  (called to find the primary)
        indirectly via ``_compute_maps``
        or directly when a primary is known. The next step, does not call it but calls
        ``_expand_hit_atom_map_by_overlap``.

        :param template:
        :param broad:
        :return:
        """
        if broad:
            self.journal.debug('Merge ligands: False. Broad search: True')
            pair_atom_maps, _ = self.get_mcs_mappings(template, self.initial_mol)
            return pair_atom_maps
        else:
            self.journal.debug('Merge ligands: False. Broad search: False')
            pair_atom_maps_t: List[IndexMap] = self._get_atom_maps(followup=self.initial_mol, hit=template,
                                                                   **self.strict_matching_mode)
            pair_atom_maps: List[Dict[int, int]] = [dict(p) for p in pair_atom_maps_t]
            return pair_atom_maps

    def _remove_negatives(self, map: Any) -> Any:
        """
        This is a hack to remove negative values from the map.
        It like flip_mapping is a tangle of options. The latter has been typing.overload'ed

        :param map:
        :return:
        """
        if len(map) == 0:  # nothing to do
            return map
        # list of something
        if isinstance(map, list) and isinstance(map[0], tuple) and isinstance(map[0][0], tuple):
            return [(h, f) for h, f in map if f >= 0 and h >= 0]
        elif isinstance(map, list) and isinstance(map[0], list):
            return [self._remove_negatives(m) for m in map]
        elif isinstance(map, list) and isinstance(map[0], dict):
            return [self._remove_negatives(m) for m in map]
        elif isinstance(map, list):
            raise ValueError('Unsupported sub map type: {}'.format(type(map[0])))
        # dict of something
        if not isinstance(map, dict):
            raise ValueError('Unsupported map type: {}'.format(type(map)))
        first_val = list(map.values())[0]
        if isinstance(map, dict) and isinstance(first_val, int):
            return {h: f for h, f in map.items() if f >= 0 and h >= 0}
        elif isinstance(map, dict) and isinstance(first_val, Iterable):
            return {k: self._remove_negatives(v) for k, v in map.items()}
        else:
            raise ValueError('Unsupported sub map type: {}'.format(type(first_val)))
