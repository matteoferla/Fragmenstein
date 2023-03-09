import json
from typing import List
from rdkit import Chem
from ..victor import Victor
from ..error import DistanceError
from rdkit.Chem import AllChem

class FreeVictor(Victor):
    """
    A test Victor: merge/place normally but no constraints are used.
    This has some weird behavior with Igor.
    """
    def post_params_step(self):
        for hydrogen in self.monster.positioned_mol.GetAtomsMatchingQuery(AllChem.AtomNumEqualsQueryAtom(1)):
            hydrogen.SetBoolProp('_Novel', True)

    def make_coordinate_constraints_for_placement(self, *args, **kwargs) -> str:
        return ''

    def make_coordinate_constraints_for_combination(self, *args, **kwargs) -> str:
        return ''


    def post_monster_step(self):
        """
        Damnatio memoriae!
        """
        self.monster.positioned_mol.ClearProp('_true_origins')
        self.monster.positioned_mol.ClearProp('_Origins')
        for atom in self.monster.positioned_mol.GetAtoms():
            atom.ClearProp('_ori_name')
            atom.ClearProp('_ori_i')
            atom.ClearProp('_x')
            atom.ClearProp('_y')
            atom.ClearProp('_z')


import json
from typing import List
from rdkit import Chem


class SingleVictor(Victor):
    """
    A test Victor: merge/place normally but the only use the hit number ``chosen``
    for the constraints.
    """
    chosen = 0

    def _calculate_combination_chem(self):
        """
        The rdkit part. Monster is used to combine the hits.
        """
        attachment = self._get_attachment_from_pdbblock() if self.is_covalent else None
        self._harmonize_warhead_combine()
        # TODO Does combine not need attachment??
        self.monster.modifications = self.modifications
        self.monster.combine(keep_all=self.monster_throw_on_discard,
                             collapse_rings=True,
                             joining_cutoff=self.joining_cutoff  # Å
                             )
        self.post_monster_step()  # empty overridable
        if self.monster_throw_on_discard and len(self.monster.unmatched):
            raise DistanceError(hits=self.monster.unmatched, distance=self.joining_cutoff)
        self.mol = self.monster.positioned_mol
        self.smiles = Chem.MolToSmiles(self.mol)
        # making folder.
        self.make_output_folder()

    def post_monster_step(self):
        """
        this empty overridable gets called after Monster does its thing
        in ``._calculate_placement_chem``,
        so before ``make_coordinate_constraints_for_placement``
        in ``._calculate_placement_chem`` gets called.
        The origins are not provided, but retrieved via ``monster.origin_from_mol``
        """
        keeper_name: str = self.hits[self.chosen].GetProp('_Name')
        self.true_hits = self.hits
        self.hits = [self.hits[self.chosen]]
        self.monster.hits = self.hits
        origins: List[List[str]] = self.monster.origin_from_mol()
        single_origins = [
            [inspiration for inspiration in atomic if keeper_name in inspiration] if isinstance(atomic, list) else []
            for atomic in origins]
        self.monster.positioned_mol.SetProp('_true_origins', json.dumps(origins))
        self.monster.positioned_mol.SetProp('_Origins', json.dumps(single_origins))
        for atom in self.monster.positioned_mol.GetAtoms():
            if atom.HasProp('_ori_name') and atom.GetProp('_ori_name') != keeper_name:
                atom.ClearProp('_ori_name')
                atom.ClearProp('_ori_i')
                atom.ClearProp('_x')
                atom.ClearProp('_y')
                atom.ClearProp('_z')