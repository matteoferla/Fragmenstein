########################################################################################################################

__doc__ = \
    """
Combine => merge/join
This merge is used by all.

Peculiarly/for historical reasons the merging is done by fragmentation and joining the fragments 
as opposed to absorption, which is what happens in the ring expansion.

    """

########################################################################################################################

from typing import Optional, Dict, List, Union
from warnings import warn

from rdkit import Chem
from rdkit.Chem import rdmolops
from ._join_neighboring import _MonsterJoinNeigh
from .bond_provenance import BondProvenance
from .positional_mapping import GPM


class _MonsterMerge(_MonsterJoinNeigh, GPM):

    def simply_merge_hits(self,
                          hits: Optional[List[Chem.Mol]] = None,
                          linked: bool = True,
                          ) -> Chem.Mol:
        """
        Recursively stick the hits together and average the positions.
        This is the monster of automerging, full-merging mapping and partial merging mapping.
        The latter however uses `partially_blend_hits` first.
        The hits are not ring-collapsed and -expanded herein.

        :param hits: optionally give a hit list, else uses the attribute ``.hits``.
        :param linked: if true the molecules are joined, else they are placed
            in the same molecule as disconnected fragments.
        :return: the rdkit.Chem.Mol object that will fill ``.scaffold``
        """
        if hits is None:
            hits = sorted(self.hits, key=lambda h: h.GetNumAtoms(), reverse=True)
        for hit in hits:
            BondProvenance.set_all_bonds(hit, 'original')
        self.journal.debug(f"Merging: {[hit.GetProp('_Name') for hit in hits]}")
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
            if linked:
                try:
                    scaffold = self.join_neighboring_mols(scaffold, fragmentanda)
                except ConnectionError:
                    self.unmatched.append(fragmentanda.GetProp("_Name"))
                    msg = f'Hit {fragmentanda.GetProp("_Name")} has no connections! Skipping!'
                    if self.throw_on_discard:
                        raise ConnectionError(msg)
                    else:
                        warn(msg)
            else:
                new_name = self._get_combined_name(scaffold, fragmentanda)
                scaffold = rdmolops.CombineMols(scaffold, fragmentanda)
                scaffold.SetProp('_Name', new_name)
        return scaffold


    def merge_pair(self, scaffold: Chem.Mol, fragmentanda: Chem.Mol, mapping: Optional = None) -> Chem.Mol:
        """
        To specify attachments use ``.merge``.
        To understand what is going on see ``.categorize``

        :param scaffold: mol to be added to.
        :param fragmentanda: mol to be fragmented
        :param mapping: see ``get_positional_mapping``. Optional in _pre_fragment_pairs
        :return:
        """
        done_already = []
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
        new_name = self._get_combined_name(scaffold, fragmentanda)
        scaffold.SetProp('_Name', new_name)
        self.keep_copy(scaffold, 'pair_merged')
        return scaffold

    def _get_combined_name(self, first_mol: Chem.Mol, second_mol: Chem.Mol):
        get_name = lambda mol: mol.GetProp('_Name') if mol.HasProp('_Name') else 'molecule'
        name_1 = get_name(first_mol)
        name_2 = get_name(second_mol)
        return f'{name_1}-{name_2}'

    # ------ merging by fragmentation  ------------------------

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

        Calls get_positional_mapping and _categorize.

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
        categories = self._categorize(fragmentanda, uniques)
        pairs = categories['pairs']
        for p in pairs:  # pairs:Dict[List[Dict]]
            for pp in pairs[p]:
                pp['idx_F'] = pp['idx']  # less ambiguous: fragmentanda index
                pp['idx_S'] = get_key(A2B_mapping, pp['idx'])  # scaffold index
        return pairs

    def _recruit_team(self, mol: Chem.Mol, starting: int, uniques: set, team: Optional[set] = None) -> set:
        if team is None:
            team = set()
        team.add(starting)
        for atom in mol.GetAtomWithIdx(starting).GetNeighbors():
            i = atom.GetIdx()
            if i in uniques and i not in team:
                team = self._recruit_team(mol, i, uniques, team)
        return team

    def _categorize(self, mol: Chem.Mol, uniques: set) -> Dict[str, Union[set, Dict]]:
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
        # if self._debug_draw: # depracated... but this could be useful...
        #     high = list(internals) + list(attachments) + list(anchors)
        #     color = {**{i: (0, 0.8, 0) for i in internals},
        #              **{i: (0, 0, 0.8) for i in attachments},
        #              **{i: (0.8, 0, 0.8) for i in anchors}}
        #     print('Purple: anchor atoms, Blue: attachments, Green: internals')
        #     self.draw_nicely(mol, highlightAtoms=high, highlightAtomColors=color)
        #     print({atom.GetIdx(): atom.GetProp('_Category') for atom in mol.GetAtoms()})
        return dict(uniques=uniques,
                    internals=internals,
                    attachments=attachments,
                    pairs=pairs,
                    dummies=dummies
                    )

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
        f = Chem.FragmentOnBonds(fragmentanda,
                                 bonds_to_frag,
                                 addDummies=False)
        frag_split = []
        fragmols = Chem.GetMolFrags(f, asMols=True, fragsMolAtomMapping=frag_split, sanitizeFrags=False)
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
        combo = Chem.RWMol(rdmolops.CombineMols(scaffold, frag))
        scaffold_anchor_index = frag_anchor_index + scaffold.GetNumAtoms()
        for detail in attachment_details:
            # scaffold_anchor_index : atom index in scaffold that needs to be added to scaffold_attachment_index
            # but was originally attached to attachment_index in fragmentanda.
            # the latter is not kept.
            attachment_index = detail['idx_F']  # fragmentanda attachment_index
            scaffold_attachment_index = detail['idx_S']  # scaffold attachment index
            bond_type = detail['type']
            combo.AddBond(scaffold_anchor_index, scaffold_attachment_index, bond_type)
            new_bond = combo.GetBondBetweenAtoms(scaffold_anchor_index, scaffold_attachment_index)
            # BondProvenance.set_bond(new_bond, '???')
            # self.transfer_ring_data(fragmentanda.GetAtomWithIdx(attachment_index),
            #                         combo.GetAtomWithIdx(scaffold_anchor_index))
        for oi, oad in zip(other_attachments, other_attachment_details):
            bond_type = oad[0]['type']
            scaffold_attachment_index = oad[0]['idx_S']
            scaffold_anchor_index = indices.index(oi) + scaffold.GetNumAtoms()
            combo.AddBond(scaffold_anchor_index, scaffold_attachment_index, bond_type)
            new_bond = combo.GetBondBetweenAtoms(scaffold_anchor_index, scaffold_attachment_index)
            # BondProvenance.set_bond(new_bond, '???')
        Chem.SanitizeMol(combo,
                         sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS +
                                     Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY,
                         catchErrors=True)
        self._prevent_two_bonds_on_dummy(combo)
        scaffold = combo.GetMol()
        return scaffold

    def _prevent_two_bonds_on_dummy(self, mol: Chem.RWMol):
        """
        The case '*(C)C' is seen legitimately in some warheads... but in most cases these are not.

        :param mol:
        :return:
        """
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != '*':
                pass
            elif len(atom.GetNeighbors()) <= 1:
                pass
            elif len(atom.GetNeighbors()) >= 2:
                self.journal.info(f'Dummy atom (idx={atom.GetIdx()}) has {len(atom.GetNeighbors())} bonds!')
                neighs = atom.GetNeighbors()
                first = neighs[0]
                for second in neighs[1:]:
                    rejected = second.GetIdx()  # that will be absorbed (deleted)
                    keeper = first.GetIdx()  # that absorbs (kept)
                    self._copy_bonding(mol, keeper, rejected)
                    self._mark_for_deletion(mol, rejected)
                self._delete_marked(mol)
                return self._prevent_two_bonds_on_dummy(mol)
