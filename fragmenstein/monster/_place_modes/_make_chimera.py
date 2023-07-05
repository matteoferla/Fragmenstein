from typing import Optional, Dict, List, Tuple, Set, Unpack, Union  # noqa cf. .legacy monkeypatch
from warnings import warn

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

from ._mappings import _MonsterMap

class _MonsterChimera(_MonsterMap):

    def make_chimera(self, template: Chem.Mol, min_mode_index=0) -> Chem.Mol:
        """
        This is to avoid extreme corner corner cases.
        E.g. here the MCS is ringMatchesRingOnly=True and AtomCompare.CompareAny,
        while for the positioning this is not the case.

        Called by full and partial blending modes.

        :return:
        """
        # get the matches
        atom_map, mode = self.get_mcs_mapping(template, self.initial_mol, min_mode_index=min_mode_index)
        follow = {**{k: str(v) for k, v in mode.items()}, 'N_atoms': len(atom_map)}
        self.journal.debug(f"scaffold-followup: {follow}")
        # make the scaffold more like the followup to avoid weird matches.
        chimera = Chem.RWMol(template)
        for scaff_ai, follow_ai in atom_map.items():
            if template.GetAtomWithIdx(scaff_ai).GetSymbol() != self.initial_mol.GetAtomWithIdx(
                    follow_ai).GetSymbol():
                v = {'F': 1, 'Br': 1, 'Cl': 1, 'H': 1, 'B': 3, 'C': 4, 'N': 3, 'O': 2, 'S': 2, 'Se': 2, 'P': 6}
                wanted = self.initial_mol.GetAtomWithIdx(follow_ai)
                if wanted.GetSymbol() == '*':  # all good then!
                    continue
                owned = template.GetAtomWithIdx(scaff_ai)
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
            chimera.UpdatePropertyCache()  # noqa
        except Chem.AtomValenceException as err:
            warn('Valance issue' + str(err))
        return chimera

    def place_from_map(self,
                       target_mol: Chem.Mol,
                       template_mol: Chem.Mol,
                       atom_map: Optional[Dict] = None,
                       random_seed=None) -> Chem.Mol:
        """
        This method places the atoms with known mapping
        and places the 'uniques' (novel) via an aligned mol (the 'sextant')
        This sextant business is a workaround for the fact that only minimised molecules can use the partial
        embedding function of RDKit.

        The template molecule may be actually two or more fragments,
        as happens for the no blending mode.
        In RDKit, the fragments within a "molecule" are not connected, but have sequential atom indices.

        :param target_mol: target mol
        :param template_mol: the template/scaffold to place the mol
        :param atom_map: something that get_mcs_mapping would return.
        :return:
        """
        # Note none of this malarkey: AllChem.MMFFOptimizeMolecule(ref)
        # prealignment
        if target_mol is None:
            target_mol = self.initial_mol
        sextant = Chem.Mol(target_mol)
        # the target (initial) mol ought to be already sanitised.
        # TODO Check if this is necessary
        # Chem.SanitizeMol(sextant)
        if random_seed:
            kwargs = dict(randomSeed=random_seed)
        else:
            kwargs = {}
        AllChem.EmbedMolecule(sextant, **kwargs)
        AllChem.MMFFOptimizeMolecule(sextant)
        ######################################################
        # mapping retrieval and sextant alignment
        # variables: atom_map sextant -> uniques
        if atom_map is None:
            atom_map, mode = self.get_mcs_mapping(target_mol, template_mol)
            msg = {**{k: str(v) for k, v in mode.items()}, 'N_atoms': len(atom_map)}
            self.journal.debug(f"followup-chimera' = {msg}")
        else:
            self.journal.debug(f'Place from map: using provided atom_map {atom_map}')
        self._add_atom_map_asProp(template_mol, atom_map)
        rdMolAlign.AlignMol(sextant, template_mol, atomMap=list(atom_map.items()), maxIters=500)
        # place atoms that have a known location
        putty = Chem.Mol(sextant)  # this is the molecule returned
        pconf = putty.GetConformer()
        chimera_conf = template_mol.GetConformer()
        uniques = set()  # unique atoms in followup
        for i in range(putty.GetNumAtoms()):
            p_atom = putty.GetAtomWithIdx(i)
            p_atom.SetDoubleProp('_Stdev', 0.)
            p_atom.SetProp('_Origin', 'none')
            if i in atom_map:
                ci = atom_map[i]
                c_atom = template_mol.GetAtomWithIdx(ci)
                if c_atom.HasProp('_Stdev'):
                    stdev = c_atom.GetDoubleProp('_Stdev')
                    origin = c_atom.GetProp('_Origin')
                    p_atom.SetDoubleProp('_Stdev', stdev)
                    p_atom.SetProp('_Origin', origin)
                pconf.SetAtomPosition(i, chimera_conf.GetAtomPosition(ci))
            else:
                uniques.add(i)
        ######################################################
        # I be using a sextant for dead reckoning!
        # variables: sextant unique team
        categories = self._categorize(sextant, uniques)
        done_already = []  # multi-attachment issue.
        for unique_idx in categories['pairs']:  # attachment unique indices
            # check the index was not done already (by virtue of a second attachment)
            if unique_idx in done_already:
                continue
            # get other attachments if any.
            team = self._recruit_team(target_mol, unique_idx, categories['uniques'])
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
            if self.attachment and list(categories['dummies']) and list(categories['dummies'])[0] in team:
                r = list(categories['dummies'])[0]
                pconf.SetAtomPosition(r, self.attachment.GetConformer().GetAtomPosition(0))
                sights.add((r, r))
            rdMolAlign.AlignMol(sextant, putty, atomMap=list(sights), maxIters=500)
            sconf = sextant.GetConformer()
            self.journal.debug(f'alignment atoms for {unique_idx} ({team}): {sights}')
            # self.draw_nicely(sextant, highlightAtoms=[a for a, b in sights])
            # copy position over
            for atom_idx in team:
                pconf.SetAtomPosition(atom_idx, sconf.GetAtomPosition(atom_idx))
            # the ring problem does not apply here but would result in rejiggling atoms.

            for other in other_attachments:
                done_already.append(other)
        # complete
        # Why is this performed: TODO Check if this is necessary
        # AllChem.SanitizeMol(putty)
        return putty  # positioned_mol

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
