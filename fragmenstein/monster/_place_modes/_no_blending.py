import operator
from typing import Optional, Dict, List, Tuple, Set, Unpack, Union  # noqa cf. .legacy monkeypatch

from rdkit import Chem
from rdkit.Chem import AllChem

from ._refine import _MonsterRefine
from ..mcs_mapping import IndexMap, flip_mapping
from ..unmerge_mapper import Unmerge


class _MonsterNone(_MonsterRefine):

    def no_blending(self, broad=False) -> None:
        """
        no merging is done. The hits are mapped individually. Not great for small fragments.
        """
        maps: Dict[str, List[Dict[int, int]]] = {}
        for template in self.hits:
            if broad:
                self.journal.debug('Merge ligands: False. Broad search: True')
                pair_atom_maps, _ = self.get_mcs_mappings(template, self.initial_mol)
                maps[template.GetProp('_Name')] = pair_atom_maps
            else:
                self.journal.debug('Merge ligands: False. Broad search: False')
                pair_atom_maps_t: List[IndexMap] = self._get_atom_maps(followup=self.initial_mol, hit=template,
                                                                       **self.strict_matching_mode)
                pair_atom_maps: List[Dict[int, int]] = [dict(p) for p in pair_atom_maps_t]
                maps[template.GetProp('_Name')] = pair_atom_maps
        # todo flip round Unmerge. Unmerge wants name to dict of followup to hit... which is backwards.
        um = Unmerge(followup=self.initial_mol,
                     mols=self.hits,
                     maps={name: list(map(flip_mapping, hit_mappings)) for name, hit_mappings in maps.items()},
                     no_discard=self.throw_on_discard)

        self.unmatched = [m.GetProp('_Name') for m in um.disregarded]
        if self.throw_on_discard and len(self.unmatched):
            raise ConnectionError(f'{self.unmatched} was rejected.')
        self.journal.debug(f'followup to scaffold {um.combined_map}')
        # ------------------ places the atoms with known mapping ------------------
        placed = self.place_from_map(target_mol=self.initial_mol,
                                     template_mol=um.combined_bonded,
                                     atom_map=um.combined_map,
                                     random_seed=self.random_seed)

        self.keep_copy(um.combined, 'scaffold')
        self.keep_copy(um.combined_bonded, 'chimera')

        alts = zip(um.combined_bonded_alternatives, um.combined_map_alternatives)
        placed_options = [self.place_from_map(target_mol=self.initial_mol,
                                              template_mol=mol,
                                              atom_map=mappa,
                                              random_seed=self.random_seed) for mol, mappa in alts]
        # ------------------ Averages the overlapping atoms ------------------
        if um.pick != -1:  # override present!
            self.journal.debug(f'override present: {um.pick}')
            self.positioned_mol = self.posthoc_refine(placed)
            self.mol_options = [self.posthoc_refine(mol) for mol in placed_options]
        else:
            mols: List[Chem.Mol] = [self.posthoc_refine(mol) for mol in [placed] + placed_options]
            self.positioned_mol: Chem.Mol = self.get_best_scoring(mols)
            mols.remove(self.positioned_mol)
            self.mol_options: List[Chem.Mol] = mols

    @classmethod
    def get_best_scoring(cls, mols: List[Chem.RWMol]) -> Chem.Mol:
        """
        Sorts molecules by how well they score w/ Merch FF
        """
        if len(mols) == 0:
            raise ValueError(f'No molecules')
        elif len(mols) == 1:
            return mols[0]
        # This is not expected to happen, but just in case
        mols = [m for m in mols if m is not None]
        scores = *map(cls.score_mol, mols),
        cls.journal.debug(f'`.get_best_scoring (unmerge)` Scores: {scores}')
        # proof that the mol has/lacks origins data:
        # for mol in mols:
        #     print('DEBUG OF LAST RESORT', mol)
        #     print(cls.origin_from_mol(cls, mol.GetMol())) # called from instance
        mol_scores = sorted(list(zip(mols, scores)), key=operator.itemgetter(1))
        return mol_scores[0][0]

    @staticmethod
    def score_mol(mol: Chem.Mol) -> float:
        """
        Scores a mol without minimising
        """
        if isinstance(mol, Chem.RWMol):
            mol = mol.GetMol()
        else:
            mol = Chem.Mol(mol)
        mol.UpdatePropertyCache()  # noqa
        Chem.SanitizeMol(mol)
        p = AllChem.MMFFGetMoleculeProperties(mol, 'MMFF94')
        if p is None:
            return float('nan')
        ff = AllChem.MMFFGetMoleculeForceField(mol, p)
        return ff.CalcEnergy()
