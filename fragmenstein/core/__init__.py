########################################################################################################################

__doc__ = \
    """
This is Fragmenstein proper. and contains the class ``Fragmenstein``.
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "0.5"
__citation__ = ""

########################################################################################################################

from typing import Dict, Union, List, Optional, Tuple
from warnings import warn
import json

import numpy as np
from collections import defaultdict, Counter

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign, rdmolops
from rdkit.Geometry.rdGeometry import Point3D

from ._utility_mixin import _FragmensteinUtil
from ._join_neighboring import _FragmensteinJoinNeighMixin
from ._collapse_ring import _FragmensteinRing
from .positional_mapping import GPM
from .unmerge_mapper import Unmerge
from .bond_provenance import BondProvenance
import itertools

##################################################################

class Fragmenstein(_FragmensteinUtil, _FragmensteinRing, GPM, _FragmensteinJoinNeighMixin):  # Unmerge is called. Not inherited.
    """
    Given a RDKit molecule and a series of hits it makes a spatially stitched together version of the initial molecule based on the hits.
    The reason is to do place the followup compound to the hits as faithfully as possible regardless of the screaming forcefields.

    * ``.scaffold`` is the combined version of the hits (rdkit.Chem.Mol object).
    * ``.scaffold_options`` are the possible scaffolds to use.
    * ``.chimera`` is the combined version of the hits, but with differing atoms made to match the followup (rdkit.Chem.Mol object).
    * ``.positioned_mol`` is the desired output (rdkit.Chem.Mol object)

    Note, the hits have to be spatially aligned —i.e. extracted from crystal structures in bond form (see. `extract_mol`).

    ``.get_positional_mapping``, which works also as a class method, creates a dictionary of mol_A atom index to mol_B atom index
    based on distance (cutoff 2&Aring;) and not MCS.

    The code works in two broad steps, first a scaffold is made, which is the combination of the hits (by position).
    Then the followup is placed. It is not embedded with constraint embed, which requires the reference molecule to have a valid geometry.
    ``.scaffold`` and ``.chimera`` and ``.positioned_mol`` absolutely do not have this.
    Novel side chains are added by aligning an optimised conformer against the closest 3-4 reference atoms.
    Note that ``.initial_mol`` is not touched. ``.positioned_mol`` may have lost some custom properties, but the atom idices are the same.

    If an atom in a Chem.Mol object is provided via ``attachment`` argument and the molecule contains a dummy atom as
    defined in the ``dummy`` class variable. Namely element R in mol file or * in string is the default.
    """

    def __init__(self, mol: Chem.Mol, hits: List[Chem.Mol], attachment: Optional[Chem.Mol] = None,
                 debug_draw: bool = False, merging_mode='partial', average_position=False):
        """
        Merging mode controls what algorithm to use.

        * 'full': a single scaffold is made
        * 'partial': multiple possible scaffolds and best is chosen
        * 'none': no merging is done. The hits are mapped individually. Not great for small fragments.
        * 'off': do nothing.

        :param mol:
        :param hits:
        :param attachment:
        :param debug_draw:
        :param merging_mode: full | partial | none | off
        """
        # starting attributes
        super().__init__()
        self.initial_mol = mol  # untouched.
        if self.initial_mol.HasSubstructMatch(self.dummy) and attachment:
            self.attachement = attachment
        elif self.initial_mol.HasSubstructMatch(self.dummy):
            warn('No attachment atom provided but dummy atom present --- ignoring.')
            self.attachement = None
        elif attachment:
            warn('Attachment atom provided but dummy atom not present --- ignoring.')
            self.attachement = None
        else:
            self.attachement = None
        # Chem.RemoveHs(self.initial_mol)
        self.hits = self.fix_hits(hits)  # list of hits
        self._debug_draw = debug_draw  # Jupyter notebook only.
        self.unmatched = []
        self.average_position = average_position
        # derived attributes
        self.scaffold = None  #: template which may have wrong elements
        self.scaffold_options = []  #: partial combined templates (merging_mode: partial)
        self.chimera = None  #: merger of hits but with atoms made to match the to-be-aligned mol
        self.positioned_mol = None  #: final molecule
        # do calculations
        if merging_mode == 'off':
            pass
        elif merging_mode == 'full':
            self.full_merging()
        elif merging_mode == 'partial':
            self.partial_merging()
        elif merging_mode == 'none_permissive' or merging_mode == 'permissive_none':
            self.no_merging(broad=True)
        elif merging_mode == 'none':
            self.no_merging()
        else:
            raise ValueError(
                f"Merging mode can only be 'full' | 'partial' | 'none' | 'none_permissive' | 'off', not '{merging_mode}'")

    def full_merging(self) -> None:
        """
        a single scaffold is made (except for ``.unmatched``)
        """
        self.scaffold_options = [self.merge_hits()]
        self.scaffold = self.posthoc_refine(self.scaffold_options[0])
        self.chimera = self.make_chimera()
        self.positioned_mol = self.place_followup()

    def partial_merging(self) -> None:
        """
        multiple possible scaffolds and best is chosen
        """
        self.scaffold_options = self.combine_hits()  # merger of hits
        unrefined_scaffold, mode_index = self.pick_best()
        used = self.scaffold.GetProp('_Name').split('-')
        self.unmatched = [h.GetProp('_Name') for h in self.hits if h.GetProp('_Name') not in used]
        self.scaffold = self.posthoc_refine(unrefined_scaffold)
        self.chimera = self.make_chimera(mode_index)
        self.positioned_mol = self.place_followup()

    def no_merging(self, broad=False) -> None:
        """
        no merging is done. The hits are mapped individually. Not great for small fragments.
        """
        maps = {}
        for template in self.hits:
            if broad:
                pair_atom_maps, _ = self.get_mcs_mappings(self.initial_mol, template)
                maps[template.GetProp('_Name')] = pair_atom_maps
            else:
                pair_atom_maps_t = self._get_atom_maps(self.initial_mol, template,
                                                       atomCompare=rdFMCS.AtomCompare.CompareElements,
                                                       bondCompare=rdFMCS.BondCompare.CompareOrder,
                                                       ringMatchesRingOnly=True,
                                                       ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
                                                       matchChiralTag=True)
                pair_atom_maps = [dict(p) for p in pair_atom_maps_t]
                maps[template.GetProp('_Name')] = pair_atom_maps
        um = Unmerge(followup=self.initial_mol,
                     mols=self.hits,
                     maps=maps,
                     no_discard=self.throw_on_discard,
                     _debug_draw=self._debug_draw)
        self.scaffold = um.combined
        full_atom_map = um.combined_map
        self.unmatched = [m.GetProp('_Name') for m in um.disregarded]
        if self.throw_on_discard and len(self.unmatched):
            raise ConnectionError(f'{self.unmatched} was rejected.')
        self.chimera = um.combined_bonded
        if self._debug_draw:
            print('followup to scaffold', full_atom_map)
            print('followup')
            self.draw_nicely(self.initial_mol)
            print('scaffold')
            self.draw_nicely(self.scaffold)
        placed = self.place_followup(atom_map=full_atom_map)
        self.positioned_mol = self.posthoc_refine(placed)

    # ========== Merging ===============================================================================================

    def merge_pair(self, scaffold: Chem.Mol, fragmentanda: Chem.Mol, mapping: Optional = None) -> Chem.Mol:
        """
        To specify attachments use ``.merge``.
        To understand what is going on see ``.categorise``

        :param scaffold: mol to be added to.
        :param fragmentanda: mol to be fragmented
        :param mapping: see ``get_positional_mapping``. Optional in _pre_fragment_pairs
        :return:
        """
        done_already = []
        if self._debug_draw:
            print('Scaffold')
            self.draw_nicely(scaffold)
            print('To be added')
            self.draw_nicely(fragmentanda)
        fp = self._pre_fragment_pairs(scaffold, fragmentanda, mapping)
        # confusingly these are hit indexed.
        for anchor_index, attachment_details in fp.items():
            # anchor index is the fragment-to-added's internal atom that attaches
            if anchor_index in done_already:
                continue
            # fix rings.
            uniques = {atom.GetIdx() for atom in fragmentanda.GetAtoms() if
                       'overlapping' not in atom.GetProp('_Category')}
            team = self._recruit_team(fragmentanda, anchor_index, uniques)
            other_attachments = list((team & set(fp.keys())) - {anchor_index})
            other_attachment_details = []
            for other in other_attachments:
                other_attachment_details.append(fp[other])
                done_already.append(other)
            scaffold = self._merge_part(scaffold, fragmentanda,
                                        anchor_index=anchor_index,
                                        attachment_details=attachment_details,
                                        other_attachments=other_attachments,
                                        other_attachment_details=other_attachment_details)
        name_A = scaffold.GetProp('_Name')
        name_B = fragmentanda.GetProp('_Name')
        scaffold.SetProp('_Name', f'{name_A}-{name_B}')
        if self._debug_draw:
            print('Merged', scaffold.GetProp('_Name'))
            self.draw_nicely(scaffold)
        return scaffold

    # ================= Combine hits ===================================================================================

    def combine_hits(self, hits: Optional[List[Chem.Mol]] = None) -> List[Chem.Mol]:
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
        for hi, hit in enumerate(hits):
            # fallback naming.
            if not hit.HasProp('_Name') or hit.GetProp('_Name').strip() == '':
                hit.SetProp('_Name', f'hit{hi}')

        ## a dodgy hit is a hit with inconsistent mapping bwteen three.
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
        if self._debug_draw:
            print('******** These combine badly')
            print(dodgy_names)
        dodgies = [hit for hit in hits if hit.GetProp('_Name') in dodgy_names]
        mergituri = [hit for hit in hits if hit.GetProp('_Name') not in dodgy_names]
        merged = self.merge_hits(mergituri)
        dodgies += [hit for hit in hits if hit.GetProp('_Name') in self.unmatched]
        self.unmatched = []
        combined_dodgies = []
        for h1, h2 in itertools.combinations(dodgies, 2):
            h_alt = Chem.Mol(h1)
            try:
                combined_dodgies.append(self.merge_pair(h_alt, h2))
            except ConnectionError:
                pass
        combinations = [merged] + dodgies + combined_dodgies
        # propagate alternatives
        while self.propagate_alternatives(combinations) != 0:
            pass
        if self._debug_draw:
            print('alternatives propagated')
        return combinations

    def propagate_alternatives(self, fewer):
        pt = Chem.GetPeriodicTable()
        new = 0
        for template in list(fewer):
            for i, atom in enumerate(template.GetAtoms()):
                if atom.HasProp('_AltSymbol'):
                    alt = Chem.Mol(template)
                    aa = alt.GetAtomWithIdx(i)
                    aa.SetAtomicNum(pt.GetAtomicNumber(atom.GetProp('_AltSymbol')))
                    aa.ClearProp('_AltSymbol')
                    atom.ClearProp('_AltSymbol')
                    fewer.append(alt)
                    new += 1
        return new

    def pick_best(self):
        if len(self.scaffold_options) == 1:
            return self.scaffold_options[0], 0
        elif len(self.scaffold_options) == 0:
            raise ValueError('No scaffolds made?!')
        else:
            mapx = {}  #: dictionary of key mol name and value tuple of maps and mode

            def template_sorter(t: List[Chem.Mol]) -> float:
                # key for sorting. requires outer scope ``maps``.
                n_atoms = len(mapx[t.GetProp('_Name')][0])
                mode = mapx[t.GetProp('_Name')][1]
                mode_i = self.matching_modes.index(mode)
                return - n_atoms - mode_i / 10

            ## get data
            # presort as this is expensive.
            for template in self.scaffold_options:
                # _get_atom_maps returns a list of alternative mappings which are lists of template to initail mol
                atom_maps = self._get_atom_maps(template, self.initial_mol,
                                                atomCompare=rdFMCS.AtomCompare.CompareElements,
                                                bondCompare=rdFMCS.BondCompare.CompareOrder,
                                                ringMatchesRingOnly=True,
                                                ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
                                                matchChiralTag=False)
                mapx[template.GetProp('_Name')] = (atom_maps, self.matching_modes[-1])
            # search properly only top 3.
            self.scaffold_options = sorted(self.scaffold_options, key=template_sorter)
            for template in self.scaffold_options[:3]:
                atom_map, mode = self.get_mcs_mapping(template, self.initial_mol)
                # get_mcs_mapping returns a dict going from template index to initial.
                mapx[template.GetProp('_Name')] = (atom_map, mode)
                if self._debug_draw:
                    print(
                        f"With {template.GetProp('_Name')}, {len(atom_map)} atoms map using mode {self.matching_modes.index(mode)}")
            ## pick best template
            self.scaffold_options = sorted(self.scaffold_options, key=template_sorter)
            ## Check if missing atoms can be explained by a different one with no overlap
            best = self.scaffold_options[0]
            ## Fuse overlaps
            # best_map = maps[best.GetProp('_Name')][0]
            # full = set(range(self.initial_mol.GetNumAtoms()))
            # present = set(best_map.values())
            # missing = full - present
            # for other in self.scaffold_options:
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

    def merge_hits(self, hits: Optional[List[Chem.Mol]] = None) -> Chem.Mol:
        """
        Recursively stick the hits together and average the positions.
        This is the core of automerging, full-merging mapping and partial merging mapping.
        The latter however uses `combine_hits` first.
        The hits are not ring-collapsed and -expanded herein.

        :param hits: optionally give a hit list, else uses the attribute ``.hits``.
        :return: the rdkit.Chem.Mol object that will fill ``.scaffold``
        """
        if hits is None:
            hits = sorted(self.hits, key=lambda h: h.GetNumAtoms(), reverse=True)
        for hit in hits:
            BondProvenance.set_all_bonds(hit, 'original')
        if self._debug_draw:
            print('Merging: ', [hit.GetProp('_Name') for hit in hits])
        scaffold = Chem.Mol(hits[0])
        # first try
        save_for_later = []
        for fragmentanda in hits[1:]:
            try:
                scaffold = self.merge_pair(scaffold, fragmentanda)
            except ConnectionError:
                save_for_later.append(fragmentanda)
        # second try
        join_later = []
        for fragmentanda in save_for_later:
            try:
                scaffold = self.merge_pair(scaffold, fragmentanda)
            except ConnectionError:
                join_later.append(fragmentanda)
        # join (last ditch)
        for fragmentanda in join_later:
            try:
                scaffold = self.join_neighboring_mols(scaffold, fragmentanda)
            except ConnectionError:
                self.unmatched.append(fragmentanda.GetProp("_Name"))
                msg = f'Hit {fragmentanda.GetProp("_Name")} has no connections! Skipping!'
                if self.throw_on_discard:
                    raise ConnectionError(msg)
                else:
                    warn(msg)
        return scaffold

    def _prevent_two_bonds_on_dummy(self, mol: Chem.RWMol):
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != '*':
                pass
            elif len(atom.GetNeighbors()) <= 1:
                pass
            elif len(atom.GetNeighbors()) >= 2:
                neighs = atom.GetNeighbors()
                for second in neighs[1:]:
                    self._absorb(mol, atom.GetIdx(), second.GetIdx())
                    mol.RemoveAtom(second.GetIdx())
                    self._prevent_two_bonds_on_dummy(mol)
                    break

    # ================= Chimera ========================================================================================

    def make_chimera(self, min_mode_index=0) -> Chem.Mol:
        """
        This is to avoid extreme corner corner cases.
        E.g. here the MCS is ringMatchesRingOnly=True and AtomCompare.CompareAny,
        while for the positioning this is not the case.

        :return:
        """
        # get the matches
        atom_map, mode = self.get_mcs_mapping(self.scaffold, self.initial_mol, min_mode_index=min_mode_index)
        self.journal.debug(f"scaffold-followup: {{**{k: str(v) for k, v in mode.items()}, 'N_atoms': len(atom_map)}}")
        if self._debug_draw:
            self.draw_nicely(self.initial_mol, highlightAtoms=atom_map.values())
        ## make the scaffold more like the followup to avoid weird matches.
        chimera = Chem.RWMol(self.scaffold)
        for scaff_ai, follow_ai in atom_map.items():
            if self.scaffold.GetAtomWithIdx(scaff_ai).GetSymbol() != self.initial_mol.GetAtomWithIdx(
                    follow_ai).GetSymbol():
                v = {'F': 1, 'Br': 1, 'Cl': 1, 'H': 1, 'B': 3, 'C': 4, 'N': 3, 'O': 2, 'S': 2, 'Se': 2, 'P': 6}
                wanted = self.initial_mol.GetAtomWithIdx(follow_ai)
                if wanted.GetSymbol() == '*':  # all good then!
                    continue
                owned = self.scaffold.GetAtomWithIdx(scaff_ai)
                diff_valance = owned.GetExplicitValence() - v[wanted.GetSymbol()]
                if wanted.GetSymbol() in ('F', 'Br', 'Cl', 'C', 'H') and diff_valance > 0:
                    continue  # cannot change this.
                elif owned.GetExplicitValence() > 4 and wanted.GetSymbol() not in ('P',):
                    continue
                else:
                    newatom = Chem.Atom(wanted)
                    stdev = chimera.GetAtomWithIdx(scaff_ai).GetDoubleProp('_Stdev')
                    newatom.SetDoubleProp('_Stdev', stdev)
                    origin = chimera.GetAtomWithIdx(scaff_ai).GetProp('_Origin')
                    newatom.SetProp('_Origin', origin)
                    chimera.ReplaceAtom(scaff_ai, newatom)
                    if diff_valance > 0:
                        chimera.GetAtomWithIdx(scaff_ai).SetFormalCharge(diff_valance)
        try:
            chimera.UpdatePropertyCache()
        except Chem.AtomValenceException as err:
            warn('Valance issue' + str(err))
        return chimera

    def place_followup(self, mol: Chem.Mol = None, atom_map: Optional[Dict] = None) -> Chem.Mol:
        """
        This method places the atoms with known mapping
        and places the 'uniques' (novel) via an aligned mol (the 'sextant')

        :param mol:
        :param atom_map: something that get_mcs_mapping would return.
        :return:
        """
        # Note none of this malarkey: AllChem.MMFFOptimizeMolecule(ref)
        # prealignment
        if mol is None:
            mol = self.initial_mol
        sextant = Chem.Mol(mol)
        Chem.SanitizeMol(sextant)
        AllChem.EmbedMolecule(sextant)
        AllChem.MMFFOptimizeMolecule(sextant)
        ######################################################
        # mapping retrieval and sextant alignment
        # variables: atom_map sextant -> uniques
        if atom_map is None:
            atom_map, mode = self.get_mcs_mapping(mol, self.chimera)
            self.journal.debug(f"followup-chimera' = {{**{k: str(v) for k, v in mode.items()}, 'N_atoms': len(atom_map)}}")
        rdMolAlign.AlignMol(sextant, self.chimera, atomMap=list(atom_map.items()), maxIters=500)
        # debug print
        if self._debug_draw:
            self.draw_nicely(mol, highlightAtoms=dict(atom_map).keys())
            self.draw_nicely(self.chimera, highlightAtoms=dict(atom_map).values())
        # place atoms that have a known location
        putty = Chem.Mol(sextant)
        pconf = putty.GetConformer()
        chimera_conf = self.chimera.GetConformer()
        uniques = set()  # unique atoms in followup
        for i in range(putty.GetNumAtoms()):
            p_atom = putty.GetAtomWithIdx(i)
            p_atom.SetDoubleProp('_Stdev', 0.)
            p_atom.SetProp('_Origin', 'none')
            if i in atom_map:
                ci = atom_map[i]
                c_atom = self.chimera.GetAtomWithIdx(ci)
                if c_atom.HasProp('_Stdev'):
                    stdev = c_atom.GetDoubleProp('_Stdev')
                    origin = c_atom.GetAtomWithIdx(ci).GetProp('_Origin')
                    p_atom.SetDoubleProp('_Stdev', stdev)
                    p_atom.SetProp('_Origin', origin)
                pconf.SetAtomPosition(i, chimera_conf.GetAtomPosition(ci))
            else:
                uniques.add(i)
        ######################################################
        # I be using a sextant for dead reckoning!
        # variables: sextant unique team
        categories = self._categorise(sextant, uniques)
        # debug print
        if self._debug_draw:
            print('internal', categories['internals'])
        done_already = []  # multi-attachment issue.
        for unique_idx in categories['pairs']:  # attachment unique indices
            # check the index was not done already (by virtue of a second attachment)
            if unique_idx in done_already:
                continue
            # get other attachments if any.
            team = self._recruit_team(mol, unique_idx, categories['uniques'])
            other_attachments = (team & set(categories['pairs'].keys())) - {unique_idx}
            sights = set()  # atoms to align against
            for att_idx in [unique_idx] + list(other_attachments):
                for pd in categories['pairs'][att_idx]:
                    first_sight = pd['idx']
                    sights.add((first_sight, first_sight))
                    neighs = [i.GetIdx() for i in sextant.GetAtomWithIdx(first_sight).GetNeighbors() if
                              i.GetIdx() not in uniques]
                    for n in neighs:
                        sights.add((n, n))
            if self.attachement and list(categories['dummies']) and list(categories['dummies'])[0] in team:
                r = list(categories['dummies'])[0]
                pconf.SetAtomPosition(r, self.attachement.GetConformer().GetAtomPosition(0))
                sights.add((r, r))
            rdMolAlign.AlignMol(sextant, putty, atomMap=list(sights), maxIters=500)
            sconf = sextant.GetConformer()
            # debug print
            if self._debug_draw:
                print(f'alignment atoms for {unique_idx} ({team}): {sights}')
                self.draw_nicely(sextant, highlightAtoms=[a for a, b in sights])
            # copy position over
            for atom_idx in team:
                pconf.SetAtomPosition(atom_idx, sconf.GetAtomPosition(atom_idx))
            # the ring problem does not apply here but would result in rejiggling atoms.

            for other in other_attachments:
                done_already.append(other)
        # complete
        AllChem.SanitizeMol(putty)
        return putty  # positioned_mol

    def _merge_part(self, scaffold: Chem.Mol, fragmentanda: Chem.Mol, anchor_index: int,
                    attachment_details: List[Dict],
                    other_attachments: List[int],
                    other_attachment_details: List[List[Dict]]) -> Chem.Mol:
        """
        This does the messy work for merge_pair.

        :param scaffold: the Chem.Mol molecule onto whose copy the fragmentanda Chem.Mol gets added
        :param fragmentanda: The other Chem.Mol molecule
        :param anchor_index: the fragment-to-added's internal atom that attaches (hit indexed)
        :param attachment_details: see `_pre_fragment_pairs` or example below fo an entry
        :type attachment_details: List[Dict]
        :param other_attachments:
        :param other_attachment_details:
        :return: a new Chem.Mol molecule

        Details object example:

            [{'idx': 5,
              'type': rdkit.Chem.rdchem.BondType.SINGLE,
              'idx_F': 5, # fragmentanda index
              'idx_S': 1  # scaffold index
              }], ...}
        """
        # get bit to add.
        bonds_to_frag = []
        for detail in attachment_details:
            attachment_index = detail['idx_F']  # fragmentanda attachment_index
            bonds_to_frag += [fragmentanda.GetBondBetweenAtoms(anchor_index, attachment_index).GetIdx()]
        bonds_to_frag += [fragmentanda.GetBondBetweenAtoms(oi, oad[0]['idx_F']).GetIdx() for oi, oad in
                          zip(other_attachments, other_attachment_details)]
        if self._debug_draw and other_attachments:
            print('ring!', other_attachments)
            print('ring!', other_attachment_details)
        f = Chem.FragmentOnBonds(fragmentanda,
                                 bonds_to_frag,
                                 addDummies=False)
        frag_split = []
        fragmols = Chem.GetMolFrags(f, asMols=True, fragsMolAtomMapping=frag_split, sanitizeFrags=False)
        if self._debug_draw:
            print('Fragment splits')
            print(frag_split)
        # Get the fragment of interest.
        ii = 0
        for mol_N, indices in enumerate(frag_split):
            if anchor_index in indices:
                break
            ii += len(indices)
        else:
            raise Exception
        frag = fragmols[mol_N]
        frag_anchor_index = indices.index(anchor_index)
        # pre-emptively fix atom ori_i
        # offset collapsed to avoid clashes.
        self.offset(frag)
        # Experimental code.
        # TODO: finish!
        # frag_atom = frag.GetAtomWithIdx(frag_anchor_index)
        # old2future = {atom.GetIntProp('_ori_i'): atom.GetIdx() + scaffold.GetNumAtoms() for atom in frag.GetAtoms()}
        # del old2future[-1] # does nothing but nice to double tap
        # if frag_atom.GetIntProp('_ori_i') == -1: #damn.
        #     for absent in self._get_mystery_ori_i(frag):
        #         old2future[absent] = scaffold_attachment_index
        # self._renumber_original_indices(frag, old2future)
        if self._debug_draw:
            print('Fragment to add')
            self.draw_nicely(frag)
        combo = Chem.RWMol(rdmolops.CombineMols(scaffold, frag))
        scaffold_anchor_index = frag_anchor_index + scaffold.GetNumAtoms()
        if self._debug_draw:
            print('Pre-merger')
            print(scaffold_anchor_index, attachment_details, anchor_index, scaffold.GetNumAtoms())
            self.draw_nicely(combo)
        for detail in attachment_details:
            # scaffold_anchor_index : atom index in scaffold that needs to be added to scaffold_attachment_index
            # but was originally attached to attachment_index in fragmentanda.
            # the latter is not kept.
            attachment_index = detail['idx_F']  # fragmentanda attachment_index
            scaffold_attachment_index = detail['idx_S']  # scaffold attachment index
            bond_type = detail['type']
            combo.AddBond(scaffold_anchor_index, scaffold_attachment_index, bond_type)
            new_bond = combo.GetBondBetweenAtoms(scaffold_anchor_index, scaffold_attachment_index)
            #BondProvenance.set_bond(new_bond, '???')
            # self.transfer_ring_data(fragmentanda.GetAtomWithIdx(attachment_index),
            #                         combo.GetAtomWithIdx(scaffold_anchor_index))
        for oi, oad in zip(other_attachments, other_attachment_details):
            bond_type = oad[0]['type']
            scaffold_attachment_index = oad[0]['idx_S']
            scaffold_anchor_index = indices.index(oi) + scaffold.GetNumAtoms()
            combo.AddBond(scaffold_anchor_index, scaffold_attachment_index, bond_type)
            new_bond = combo.GetBondBetweenAtoms(scaffold_anchor_index, scaffold_attachment_index)
            #BondProvenance.set_bond(new_bond, '???')
            if self._debug_draw:
                print(
                    f"Added additional {bond_type.name} bond between {scaffold_attachment_index} and {scaffold_anchor_index} " + \
                    f"(formerly {indices.index(oi)})")
        Chem.SanitizeMol(combo,
                         sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS +
                                     Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY,
                         catchErrors=True)
        if self._debug_draw:
            print('Merged')
            self.draw_nicely(combo)
        self._prevent_two_bonds_on_dummy(combo)
        scaffold = combo.GetMol()
        return scaffold

    def transfer_ring_data(self, donor: Chem.Atom, acceptor: Chem.Atom):
        """
        Transfer the info if a ringcore atom.

        :param donor:
        :param acceptor:
        :return:
        """
        # if donor.GetIntProp('_ori_i') == -1:
        #     data = donor
        pass

    def _pre_fragment_pairs(self, scaffold: Chem.Mol, fragmentanda: Chem.Mol, A2B_mapping: Optional = None) \
            -> Dict[int, List[Dict]]:
        """
        Returns

            {4: [{'idx': 5,
                   'type': rdkit.Chem.rdchem.BondType.SINGLE,
                   'idx_F': 5,
                   'idx_S': 1}], ...}

        which is slight more than {5: [{'idx': 4, 'type': rdkit.Chem.rdchem.BondType.SINGLE}], ... from categories

        idx_F: fragmentanda index
        idx_S: scaffold index

        required for self.merge, the key is the index of anchoring atom.

        Calls get_positional_mapping and _categorise.

        :param scaffold: mol to be added to.
        :param fragmentanda: mol to be fragmented
        :param A2B_mapping: see ``get_positional_mapping``
        :return:
        """
        # get A2B mapping
        if A2B_mapping is None:
            A2B_mapping = self.get_positional_mapping(scaffold, fragmentanda)
        get_key = lambda d, v: list(d.keys())[list(d.values()).index(v)]
        if len(A2B_mapping) == 0:
            raise ConnectionError('No overlap!')
        # store alternative atom symbols.
        for si, fi in A2B_mapping.items():
            sa = scaffold.GetAtomWithIdx(si)
            sn = sa.GetSymbol()
            fn = fragmentanda.GetAtomWithIdx(fi).GetSymbol()
            if sn != fn:
                sa.SetProp('_AltSymbol', fn)
        # prepare.
        uniques = set(range(fragmentanda.GetNumAtoms())) - set(A2B_mapping.values())
        categories = self._categorise(fragmentanda, uniques)
        pairs = categories['pairs']
        for p in pairs:  # pairs:Dict[List[Dict]]
            for pp in pairs[p]:
                pp['idx_F'] = pp['idx']  # less ambiguous: fragmentanda index
                pp['idx_S'] = get_key(A2B_mapping, pp['idx'])  # scaffold index
        return pairs

    def _categorise(self, mol: Chem.Mol, uniques: set) -> Dict[str, Union[set, Dict]]:
        """
        What do the novel atoms do in terms of connectivity.
        Complicated dict output (called ``categories`` in the methods). Really ought to be SetProp of the atoms.

        * ``uniques`` are set of atoms to classify on
        * ``internals`` are unique atoms that are connected solely to unique atoms
        * ``attachments`` are non-unique atoms to which a unique atom connects
        * ``pairs`` is a dict of unique atom idx --> dict of ``idx`` --> attachment idx and ``type`` bond type.

        :param mol: molecule to describe
        :param uniques: set of indices that are new to this molecule
        :return:
        """
        #
        pairs = {}
        internals = set()
        attachments = set()
        dummies = set()
        for i in uniques:  # novel atoms
            unique_atom = mol.GetAtomWithIdx(i)
            if unique_atom.GetSymbol() == self.dummy_symbol:
                dummies.add(i)
            neighbours = {n.GetIdx() for n in unique_atom.GetNeighbors()}
            if len(neighbours - uniques) == 0:  # unlessone of the connections is not unique.
                internals.add(i)
            else:
                i_attached = neighbours - uniques
                attachments |= i_attached
                pairs[i] = [{'idx': j,
                             'type': mol.GetBondBetweenAtoms(i, j).GetBondType()} for j in i_attached]
        anchors = uniques - internals
        # store for safekeeping
        for atom in mol.GetAtoms():
            i = atom.GetIdx()
            if i in internals:  # novel and not connected
                atom.SetProp('_Category', 'internal')
            elif i in attachments:  # not-novel but connected
                atom.SetProp('_Category', 'overlapping-attachment')
            elif i in pairs:  # dict not set tho
                atom.SetProp('_Category', 'internal-attachment')
            else:  # overlapping
                atom.SetProp('_Category', 'overlapping')
        if self._debug_draw:
            high = list(internals) + list(attachments) + list(anchors)
            color = {**{i: (0, 0.8, 0) for i in internals},
                     **{i: (0, 0, 0.8) for i in attachments},
                     **{i: (0.8, 0, 0.8) for i in anchors}}
            print('Purple: anchor atoms, Blue: attachments, Green: internals')
            self.draw_nicely(mol, highlightAtoms=high, highlightAtomColors=color)
            print({atom.GetIdx(): atom.GetProp('_Category') for atom in mol.GetAtoms()})
        return dict(uniques=uniques,
                    internals=internals,
                    attachments=attachments,
                    pairs=pairs,
                    dummies=dummies
                    )

    # ========= Other ==================================================================================================

    def fix_hits(self, hits: List[Chem.Mol]) -> List[Chem.Mol]:
        for hi, hit in enumerate(hits):
            if isinstance(hit, str):
                warn(f'Hit {hi} is a string ({hit}). This route is not the intended way. Trying to read it.')
                if '.mol' in hit or '.mdf' in hit:
                    hits[hi] = Chem.MolFromMolFile(hit)
                elif '.pdb' in hit:
                    hits[hi] = Chem.MolFromPDBFile(hit)
                else:
                    raise ValueError(f'Hit {hit} is not a Mol file.')
            elif isinstance(hit, Chem.Mol):
                pass
            else:
                raise ValueError(f'Hit has to be a Chem.Mol! not {type(hit)}')
            # fallback naming.
            if not hit.HasProp('_Name') or hit.GetProp('_Name').strip() == '':
                hit.SetProp('_Name', f'hit{hi}')
            self.store_positions(hit)
        return hits

    def posthoc_refine(self, scaffold, indices: Optional[List[int]] = None) -> Chem.Mol:
        """
        Averages the overlapping atoms.

        :param scaffold:
        :return:
        """
        if indices is None:
            indices = list(range(scaffold.GetNumAtoms()))
        refined = Chem.RWMol(scaffold)
        refconf = refined.GetConformer()
        positions = defaultdict(list)  # coordinates
        equivalence = defaultdict(list)  # atom indices of hits.
        for h in self.hits:
            if h.GetProp('_Name') in self.unmatched:
                continue
            hc = h.GetConformer()
            for k, v in self.get_positional_mapping(scaffold, h).items():
                positions[k].append([hc.GetAtomPosition(v).x, hc.GetAtomPosition(v).y, hc.GetAtomPosition(v).z])
                equivalence[k].append(f'{h.GetProp("_Name")}.{v}')
        for i in range(scaffold.GetNumAtoms()):
            if i not in indices:
                continue
            elif len(positions[i]) == 0:
                refined.GetAtomWithIdx(i).SetDoubleProp('_Stdev', 0.)
                refined.GetAtomWithIdx(i).SetDoubleProp('_Max', 0.)
                refined.GetAtomWithIdx(i).SetProp('_Origin', 'none')
                # warn(f'Atom {i}  {scaffold.GetAtomWithIdx(i).GetSymbol}/{refined.GetAtomWithIdx(i).GetSymbol} '+ \
                #     'in scaffold that has no positions.')
            else:
                p = np.mean(np.array(positions[i]), axis=0).astype(float)
                # sd = np.mean(np.std(np.array(positions[i]), axis=0)).astype(float)
                ds = [np.linalg.norm(p - pi) for pi in positions[i]]
                sd = np.std(ds)
                md = np.max(ds)
                refined.GetAtomWithIdx(i).SetProp('_Origin', json.dumps(equivalence[i]))
                refined.GetAtomWithIdx(i).SetDoubleProp('_Stdev', sd)
                refined.GetAtomWithIdx(i).SetDoubleProp('_Max', md)
                if self.average_position:
                    refconf.SetAtomPosition(i, Point3D(p[0], p[1], p[2]))
        Chem.SanitizeMol(refined,
                         sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS +
                                     Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY,
                         catchErrors=True)
        return refined

    def get_mcs_mappings(self, molA, molB, min_mode_index: int = 0) -> Tuple[List[Dict[int, int]], dict]:
        """
        This is a weird method. It does a strict MCS match.
        And then it uses laxer searches and finds the case where a lax search includes the strict search.

        :param molA: query molecule
        :param molB: target/ref molecule
        :param min_mode_index: the lowest index to try (opt. speed reasons)
        :return: mappings and mode
        """
        strict = self._get_atom_maps(molA, molB, atomCompare=rdFMCS.AtomCompare.CompareElements,
                                     bondCompare=rdFMCS.BondCompare.CompareOrder,
                                     ringMatchesRingOnly=True,
                                     ringCompare=rdFMCS.RingCompare.PermissiveRingFusion,
                                     matchChiralTag=True)
        for i, mode in enumerate(self.matching_modes):
            if i < min_mode_index:
                continue
            if self._debug_draw:
                print(f'Trying mode {mode}')
            lax = self._get_atom_maps(molA, molB, **mode)
            # remove the lax matches that disobey
            neolax = [l for l in lax if any([len(set(s) - set(l)) == 0 for s in strict])]
            if len(neolax) == 0:
                continue
            else:
                return [dict(n) for n in neolax], mode
        else:
            raise ValueError('This is chemically impossible: nothing matches in the MCS step ' +\
                             f'({len(self.matching_modes)}')

    def get_mcs_mapping(self, molA, molB, min_mode_index: int = 0) -> Tuple[Dict[int, int], dict]:
        """
        This is a weird method. It does a strict MCS match.
        And then it uses laxer searches and finds the case where a lax search includes the strict search.

        :param molA: query molecule
        :param molB: target/ref molecule
        :param min_mode_index: the lowest index to try (opt. speed reasons)
        :return: mapping and mode
        """
        ms, mode = self.get_mcs_mappings(molA, molB, min_mode_index)
        return ms[0], mode

    def _get_atom_maps(self, molA, molB, **mode) -> List[List[Tuple[int, int]]]:
        mcs = rdFMCS.FindMCS([molA, molB], **mode)
        common = Chem.MolFromSmarts(mcs.smartsString)
        matches = []
        # prevent a dummy to match a non-dummy, which can happen when the mode is super lax.
        is_dummy = lambda mol, at: mol.GetAtomWithIdx(at).GetSymbol() == '*'
        all_bar_dummy = lambda Aat, Bat: (is_dummy(molA, Aat) and is_dummy(molB, Bat)) or not (
                is_dummy(molA, Aat) or is_dummy(molB, Bat))
        for molA_match in molA.GetSubstructMatches(common, uniquify=False):
            for molB_match in molB.GetSubstructMatches(common, uniquify=False):
                matches.append([(molA_at, molB_at) for molA_at, molB_at in zip(molA_match, molB_match) if
                                all_bar_dummy(molA_at, molB_at)])
        # you can map two toluenes 4 ways, but two are repeats.
        matches = set([tuple(sorted(m, key=lambda i: i[0])) for m in matches])
        return matches

    def _get_atom_map(self, molA, molB, **mode) -> List[Tuple[int, int]]:
        return self._get_atom_maps(molA, molB, **mode)[0]

    def _recruit_team(self, mol: Chem.Mol, starting: int, uniques: set, team: Optional[set] = None) -> set:
        if team is None:
            team = set()
        team.add(starting)
        for atom in mol.GetAtomWithIdx(starting).GetNeighbors():
            i = atom.GetIdx()
            if i in uniques and i not in team:
                team = self._recruit_team(mol, i, uniques, team)
        return team

    def pretweak(self) -> None:
        """
        What if the fragments were prealigned slightly? Really bad things.

        :return:
        """
        warn('This method is unreliable. Do not use it')
        ref = self.hits[0]
        for target in self.hits[1:]:
            A2B = list(self.get_positional_mapping(target, ref, 0.5).items())
            if A2B:
                rdMolAlign.AlignMol(target, ref, atomMap=A2B, maxIters=500)
            else:
                warn(f'No overlap? {A2B}')
