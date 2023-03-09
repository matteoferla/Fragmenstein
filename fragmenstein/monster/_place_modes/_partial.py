import warnings

from rdkit.Chem import rdmolops

import itertools, operator
import json
from collections import Counter
from collections import defaultdict
from typing import Optional, Dict, List, Tuple, Set, Unpack, Union  # noqa cf. .legacy monkeypatch
from warnings import warn

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS, rdMolAlign, rdmolops
from rdkit.Chem import rdmolops
from rdkit.Geometry.rdGeometry import Point3D

from .._communal import _MonsterCommunal
from .._merge import _MonsterMerge
from ..unmerge_mapper import Unmerge
from ..mcs_mapping import SpecialCompareAtoms, IndexMap, ExtendedFMCSMode, transmute_FindMCS_parameters, flip_mapping
from ._refine import _MonsterRefine
from ...error import FragmensteinError

class _MonsterPartial(_MonsterRefine):
    # placement dependent methods

    def partial_blending(self) -> None:
        """
        multiple possible scaffolds for placement and best is chosen
        """
        self.mol_options = self.partially_blend_hits()  # merger of hits
        unrefined_scaffold, mode_index = self.pick_best()
        used = unrefined_scaffold.GetProp('_Name').split('-')
        self.unmatched = [h.GetProp('_Name') for h in self.hits if h.GetProp('_Name') not in used]
        scaffold = self.posthoc_refine(unrefined_scaffold)
        chimera = self.make_chimera(scaffold, mode_index)

        # atom_map is filled via `self.get_mcs_mapping(target_mol, template_mol)`:
        self.positioned_mol = self.place_from_map(target_mol=self.positioned_mol,
                                                  template_mol=chimera,
                                                  atom_map=None,
                                                  random_seed=self.random_seed)
        self.keep_copy(scaffold, 'scaffold')
        self.keep_copy(chimera, 'chimera')

    def partially_blend_hits(self, hits: Optional[List[Chem.Mol]] = None) -> List[Chem.Mol]:
        """
        This is the partial merge algorithm, wherein the hits are attempted to be combined.
        If the combination is bad. It will not be combined.
        Returning a list of possible options.
        These will have the atoms changed too.

        :param hits:
        :param distance:
        :return:
        """

        if hits is None:
            hits = sorted(self.hits, key=lambda h: h.GetNumAtoms(), reverse=True)
        # this predates fix_hits. Check if it is still needed.
        for hi, hit in enumerate(hits):
            # fallback naming.
            if not hit.HasProp('_Name') or hit.GetProp('_Name').strip() == '':
                hit.SetProp('_Name', f'hit{hi}')

        ## a dodgy hit is a hit with inconsistent mapping between three.
        def get_dodgies(skippers):
            dodgy = []
            for hit0, hit1, hit2 in itertools.combinations(hits, 3):
                hn0 = hit0.GetProp('_Name')
                hn1 = hit1.GetProp('_Name')
                hn2 = hit2.GetProp('_Name')
                if any([hit in skippers for hit in (hn0, hn1, hn2)]):
                    continue
                for a, b in inter_mapping[(hn0, hn1)].items():
                    if a in inter_mapping[(hn0, hn2)] and b in inter_mapping[(hn1, hn2)]:
                        if inter_mapping[(hn0, hn2)][a] != inter_mapping[(hn1, hn2)][b]:
                            # TODO: THIS IS A BAD OPTION:
                            # if all([m.GetAtomWithIdx(i).IsInRing() for m, i in ((hit0, a),
                            #                                                     (hit1, b),
                            #                                                     (hit2, inter_mapping[(hn0, hn2)][a]),
                            #                                                     (hit2, inter_mapping[(hn1, hn2)][b]))]):
                            #     pass
                            # else:
                            dodgy.extend((hn0, hn1, hn2))
            d = Counter(dodgy).most_common()
            if dodgy:
                return get_dodgies(skippers=skippers + [d[0][0]])
            else:
                return skippers

        inter_mapping = {}
        for h1, h2 in itertools.combinations(hits, 2):
            inter_mapping[(h1.GetProp('_Name'), h2.GetProp('_Name'))] = self.get_positional_mapping(h1, h2)
        dodgy_names = get_dodgies([])
        warn(f'These combiend badly: {dodgy_names}')
        dodgies = [hit for hit in hits if hit.GetProp('_Name') in dodgy_names]
        mergituri = [hit for hit in hits if hit.GetProp('_Name') not in dodgy_names]
        merged = self.simply_merge_hits(mergituri)
        dodgies += [hit for hit in hits if hit.GetProp('_Name') in self.unmatched]
        self.unmatched = []
        combined_dodgies = []
        for h1, h2 in itertools.combinations(dodgies, 2):
            h_alt = Chem.Mol(h1)
            try:
                combined_dodgies.append(self.merge_pair(h_alt, h2))
            except FragmensteinError:
                pass
        combinations = [merged] + dodgies + combined_dodgies
        # propagate alternatives
        while self.propagate_alternatives(combinations) != 0:
            pass
        return combinations

    def pick_best(self) -> Tuple[Chem.Mol, int]:
        """
        Method for partial merging for placement

        :return: unrefined_scaffold, mode_index
        """
        if len(self.mol_options) == 1:
            return self.mol_options[0], 0
        elif len(self.mol_options) == 0:
            raise ValueError('No scaffolds made?!')
        else:
            mapx = {}  #: dictionary of key mol name and value tuple of maps and mode

            def template_sorter(t: List[Chem.Mol]) -> float:
                # key for sorting. requires outer scope ``maps``.
                self.journal.debug(f'pick_best (partial merging): Sorting {len(t)} templates: {t}')
                n_atoms = len(mapx[t[0].GetProp('_Name')])
                mode = mapx[t[1].GetProp('_Name')]
                mode_i = self.matching_modes.index(mode)
                return - n_atoms - mode_i / 10

            ## get data
            # presort as this is expensive.
            for template in self.mol_options:
                # _get_atom_maps returns a list of alternative mappings which are lists of template to initail mol
                # todo flip around
                atom_maps: List[IndexMap] = self._get_atom_maps(followup=self.initial_mol, hit=template,
                                                                atomCompare=rdFMCS.AtomCompare.CompareElements,
                                                                bondCompare=rdFMCS.BondCompare.CompareOrder,
                                                                ringMatchesRingOnly=True,
                                                                ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
                                                                matchChiralTag=False)
                atom_maps: List[IndexMap] = flip_mapping(atom_maps)
                mapx[template.GetProp('_Name')] = (atom_maps, self.matching_modes[-1])
            # search properly only top 3.
            self.mol_options = sorted(self.mol_options, key=template_sorter)
            for template in self.mol_options[:3]:
                atom_map, mode = self.get_mcs_mapping(template, self.initial_mol)
                # get_mcs_mapping returns a dict going from template index to initial.
                mapx[template.GetProp('_Name')] = (atom_map, mode)
                self.journal.debug(f"With {template.GetProp('_Name')}, " + \
                                   "{len(atom_map)} atoms map using mode {self.matching_modes.index(mode)}")
            ## pick best template
            self.mol_options = sorted(self.mol_options, key=template_sorter)
            ## Check if missing atoms can be explained by a different one with no overlap
            best = self.mol_options[0]
            ## Fuse overlaps
            # best_map = maps[best.GetProp('_Name')][0]
            # full = set(range(self.initial_mol.GetNumAtoms()))
            # present = set(best_map.values())
            # missing = full - present
            # for other in self.mol_options:
            #     other_map = maps[other.GetProp('_Name')][0]
            #     found = set(other_map.values())
            #     if len(found) > 6 and len(present & found) == 0: # more than just a ring and no overlap
            #         fusion = self._fuse(best, other, best_map, other_map)
            return best, self.matching_modes.index(mapx[best.GetProp('_Name')][1])

    # def _fuse(self, mol_A: Chem.Mol, mol_B: Chem.Mol, map_A: Dict[int, int], map_B: Dict[int, int]) -> Chem.Mol:
    #     """
    #     Merge two compounds... but that are unlinked, using the followup as a guide.
    #     Conceptually different but overlapping is join_neighboring_mols
    #
    #     :param mol_A:
    #     :param mol_B:
    #     :param map_A:
    #     :param map_B:
    #     :return:
    #     """
    #     # No longer needed.
    #     fusion = Chem.RwMol(Chem.CombineMols(mol_A, mol_B))
    #     t = mol_A.GetNumAtoms()
    #     new_map_B = {k+t: v for k, v in map_B.items()}
    #     full = set(range(self.initial_mol.GetNumAtoms()))
    #     present_A = set(map_A.values())
    #     present_B = set(map_B.values())
    #
    #     def find_route(n):
    #         if n in present_A:
    #             return None
    #         elif n in present_B:
    #             return n
    #         else:
    #             path_raw = {m: find_route(m) for m in self.initial_mol.GetAtomWithIdx(n).GetNeighbors()}
    #             path = {i: path_raw[i] for i in path_raw if path_raw[i] is not None}
    #             if len(path) == 0:
    #                 return None
    #             else:
    #                 return {n: path}
