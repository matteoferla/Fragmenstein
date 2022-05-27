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
              merging_mode: str = 'expansion'):
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
        :param attachment:
        :param custom_map:
        :param merging_mode:
        :return:
        """
        if not mol.HasProp('_Name'):
            mol.SetProp('_Name', 'followup')
        self.initial_mol, self.attachment = self._parse_mol_for_place(mol, attachment)
        if custom_map:
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
        for hit in self.hits: #: Chem.Mol
            hit_name:str = hit.GetProp('_Name')
            hit_names.append(hit_name)
            if hit_name not in self.custom_map:
                self.custom_map[hit_name] = {}
        # check for incorrect keys
        for hit_name in self.custom_map:
            if hit_name not in hit_names:
                raise ValueError(f"Custom map contains key '{hit_name}' which is not in hits ({hit_names}).")

