########################################################################################################################

__doc__ = \
    """
Place followup
    """

########################################################################################################################

from typing import Optional, Dict, List
from warnings import warn

from rdkit import Chem

from ._place_modes import _MonsterBlend
from .mcs_mapping import IndexMap

class _MonsterPlace(_MonsterBlend):
    """
    * _MonsterBase
    * _MonsterCommunal
    * _MonsterTracker
    * _MonsterMerge
    * _MonsterCombine
    * _MonsterMap
    * _MonsterChimera
    * _MonsterRefine
    * (_MonsterFull, _MonsterPartial, _MonsterNone) >
    * _MonsterBlend
    * _MonsterPlace
    """

    def place(self,
              mol: Chem.Mol,
              attachment: Optional[Chem.Mol] = None,
              custom_map: Optional[Dict[str, Dict[int, int]]] = None,
              merging_mode: str = 'expansion',
              enforce_warhead_mapping:bool=True):
        """
        Positioned a given mol based on the hits. (Main entrypoint)
        accepts the argument `merging_mode`, by default it is "expansion",
        but was "permissive_none",
        which call `.by_expansion` and `.no_blending(broad=True)` respectively.
        "off" (does nothing except fill the attribute ``initial_mol``),
        "full" (`.full_blending()`),
        "partial" (`.partial_blending()`)
        and "none" (`.no_blending()`, but less thorough)
        are accepted.

        :param mol:
        :param attachment: This the SG of the cysteine if covalent
        :param custom_map: Dict of hit_name to Dict of hit_idx to followup_idx
        :param merging_mode:
        :return:
        """
        if not mol.HasProp('_Name'):
            mol.SetProp('_Name', 'followup')
        assert not any([mol.GetProp('_Name') == h.GetProp('_Name') for h in self.hits]), 'Placement has the same name as a hit!'
        self.initial_mol, self.attachment = self._parse_mol_for_place(mol, attachment)
        custom_map = self.fix_custom_map(custom_map)
        if enforce_warhead_mapping and self.attachment:
            self.journal.debug('Enforcing warhead mapping')
            custom_map: Dict[str, Dict[int, int]] = self._add_warhead_mapping(custom_map)
        self.custom_map: Dict[str, Dict[int, int]] = custom_map
        # do the mol names match?
        self._validate_custom_map()
        # Reset
        self.unmatched = []
        self.mol_options = []
        # do calculations
        self.merging_mode = merging_mode
        if merging_mode == 'off':
            pass
        elif merging_mode == 'full':
            self.full_blending()
        elif merging_mode == 'partial':
            self.partial_blending()
        elif merging_mode == 'none_permissive' or merging_mode == 'permissive_none':
            self.no_blending(broad=True)
        elif merging_mode == 'none':
            self.no_blending()
        elif merging_mode == 'expansion':
            self.by_expansion()
        else:
            valid_modes = ('full', 'partial', 'none', 'none_permissive', 'off', 'expansion')
            raise ValueError(
                f"Merging mode can only be {'| '.join(valid_modes)}, not '{merging_mode}'")
        return self

    def place_smiles(self,
                     smiles: str,
                     long_name: Optional[str]=None,
                     **kwargs):
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp('_Name', long_name if long_name else 'followup')
        self.place(mol=mol, **kwargs)
        return self

    def _parse_mol_for_place(self,
                             mol: Chem.Mol,
                             attachment: Optional[Chem.Mol] = None):
        """
        Atachment is a single atom molecule, say the SG atom of a CYS, onto which the mol is placed.
        :param mol:
        :param attachment:
        :return:
        """
        # ------------- store mol ---------------------------------------
        if mol.HasSubstructMatch(self.dummy) and attachment:
            pass
        elif mol.HasSubstructMatch(self.dummy):
            warn('No attachment atom provided but dummy atom present --- ignoring.')
            attachment = None
        elif attachment:
            warn('Attachment atom provided but dummy atom not present --- ignoring.')
            attachment = None
        else:
            attachment = None
        return mol, attachment

    def _validate_custom_map(self):
        """
        What on earth has the user submitted as custom_map?

        :return:
        """
        # add missing keys
        hit_names:List[str] = []
        for hit in self.hits:  #: Chem.Mol
            hit_name:str = hit.GetProp('_Name')
            hit_names.append(hit_name)
            if hit_name not in self.custom_map:
                self.custom_map[hit_name] = {}
        # check for incorrect keys
        for hit_name in self.custom_map:
            if hit_name not in hit_names:
                raise ValueError(f"Custom map contains key '{hit_name}' which is not in hits ({hit_names}).")

    def _add_warhead_mapping(self, custom_map: Dict[str, IndexMap]) -> Dict[str, IndexMap]:
        """
        Add the warhead mapping to the custom_map.
        Does not use the warhead definition as it may be missing.

        :param custom_map:
        :return:
        """
        q: Chem.rdchem.QueryAtom = Chem.rdqueries.AtomNumEqualsQueryAtom(0)
        dummies = list(self.initial_mol.GetAtomsMatchingQuery(q))
        if len(dummies) == 0:
            self.journal.warning('The provided molecule is not covalent...')
            return custom_map
        if len(dummies) > 1:
            self.journal.info('More than one dummy atom found: cannot enforce warhead mapping.')
            return custom_map
        followup_dummy = dummies[0]

        def get_neighs(atom: Chem.Atom, old: Chem.Atom) -> List[Chem.Atom]:
            """get hetero neighs sorted by weight, except for the previous one"""
            neighs = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1 and n.GetIdx() != old.GetIdx()]
            return sorted(neighs, key=lambda n: n.GetAtomicNum())

        for hit in self.hits:
            dummies = list(hit.GetAtomsMatchingQuery(q))
            if len(dummies) != 1:
                continue
            hit_dummy = dummies[0]
            hit_name = hit.GetProp('_Name')
            if hit_name not in custom_map:
                custom_map[hit_name] = {}
            # the hit and followup have a dummy atom:
            custom_map[hit_name][hit_dummy.GetIdx()] = followup_dummy.GetIdx()
            followup_atom = followup_dummy
            hit_atom = hit_dummy
            followup_old, hit_old = followup_atom, hit_atom
            while True:
                # crawl up the network
                followup_hetero_neighs = get_neighs(followup_atom, followup_old)
                hit_hetero_neighs = get_neighs(hit_atom, hit_old)
                if len(followup_hetero_neighs) == 1 and len(hit_hetero_neighs) == 1:
                    followup_old, hit_old = followup_atom, hit_atom
                    followup_atom = followup_hetero_neighs[0]
                    hit_atom = hit_hetero_neighs[0]
                    custom_map[hit_name][hit_atom.GetIdx()] = followup_atom.GetIdx()
                elif len(followup_hetero_neighs) == 2 and len(hit_hetero_neighs) == 2 and \
                        followup_hetero_neighs[0].GetSymbol() == 'C' and \
                        hit_hetero_neighs[0].GetSymbol() == 'C' and \
                        followup_hetero_neighs[1].GetSymbol() in ('O', 'N') and \
                        hit_hetero_neighs[1].GetSymbol() in ('O', 'N'):
                    custom_map[hit_name][hit_hetero_neighs[0].GetIdx()] = followup_hetero_neighs[0].GetIdx()
                    custom_map[hit_name][hit_hetero_neighs[1].GetIdx()] = followup_hetero_neighs[1].GetIdx()
                    followup_old, hit_old = followup_atom, hit_atom
                    followup_atom, hit_atom = followup_hetero_neighs[0], hit_hetero_neighs[0]
                else:
                    break
        return custom_map
