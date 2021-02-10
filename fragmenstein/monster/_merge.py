########################################################################################################################

__doc__ = \
    """
Combine = merge/join
    """

########################################################################################################################

from typing import Optional, Dict, List, Union
from warnings import warn

from rdkit import Chem

from ._join_neighboring import _MonsterJoinNeigh
from .bond_provenance import BondProvenance
from .positional_mapping import GPM


class _MonsterMerge(_MonsterJoinNeigh, GPM):
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
        self.keep_copy(scaffold, 'pair_merged')
        return scaffold

    def simply_merge_hits(self, hits: Optional[List[Chem.Mol]] = None) -> Chem.Mol:
        """
        Recursively stick the hits together and average the positions.
        This is the monster of automerging, full-merging mapping and partial merging mapping.
        The latter however uses `partially_blend_hits` first.
        The hits are not ring-collapsed and -expanded herein.

        :param hits: optionally give a hit list, else uses the attribute ``.hits``.
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

    def _recruit_team(self, mol: Chem.Mol, starting: int, uniques: set, team: Optional[set] = None) -> set:
        if team is None:
            team = set()
        team.add(starting)
        for atom in mol.GetAtomWithIdx(starting).GetNeighbors():
            i = atom.GetIdx()
            if i in uniques and i not in team:
                team = self._recruit_team(mol, i, uniques, team)
        return team

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