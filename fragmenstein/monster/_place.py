########################################################################################################################

__doc__ = \
    """
Place followup
    """

########################################################################################################################

from typing import Optional, Dict
from warnings import warn

from rdkit import Chem

from ._blend_place import _MonsterBlend  # inherits _MonsterCommunal <- _MonsterBase


class _MonsterPlace(_MonsterBlend):

    def place(self,
              mol: Chem.Mol,
              attachment: Optional[Chem.Mol] = None,
              custom_map: Optional[Dict[str, Dict[int, int]]] = None,
              merging_mode: str = 'none_permissive'):
        """
        Positioned a given mol based on the hits. (Main entrypoint)
        accepts the argument `merging_mode`, by default it is "permissive_none",
        which calls `.no_blending(broad=True)`,
        but "off" (does nothing except fill the attribute ``initial_mol``),
        "full" (`.full_blending()`),
        "partial" (`.partial_blending()`)
        and "none" (`.no_blending()`)
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
        else:
            valid_modes = ('full', 'partial', 'none', 'none_permissive', 'off')
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
