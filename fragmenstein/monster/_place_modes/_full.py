from typing import Optional, Dict, List, Tuple, Set, Unpack, Union  # noqa cf. .legacy monkeypatch

from typing import Optional, Dict, List, Tuple, Set, Unpack, Union  # noqa cf. .legacy monkeypatch

from ._refine import _MonsterRefine


class _MonsterFull(_MonsterRefine):
    # placement dependent methods

    # @classmethod #why was this a classmethod
    def full_blending(self) -> None:
        """
        a single scaffold is made (except for ``.unmatched``)
        """
        self.mol_options = [self.simply_merge_hits()]
        scaffold = self.posthoc_refine(self.mol_options[0])
        chimera = self.make_chimera(scaffold)

        # atom_map is filled via `self.get_mcs_mapping(target_mol, template_mol)`:
        self.positioned_mol = self.place_from_map(target_mol=self.initial_mol,
                                                  template_mol=chimera,
                                                  atom_map=None,
                                                  random_seed=self.random_seed)
        self.keep_copy(scaffold, 'scaffold')
        self.keep_copy(chimera, 'chimera')