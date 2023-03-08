import json
from typing import List
from rdkit import Chem
from ..victor import Victor

class FreeVictor(Victor):
    """
    A test Victor: merge/place normally but no constraints are used.
    """

    def make_coordinate_constraints_for_placement(self, *args, **kwargs) -> str:
        return ''

    def make_coordinate_constraints_for_combination(self, *args, **kwargs) -> str:
        return ''


class SingleVictor(Victor):
    """
    A test Victor: merge/place normally but the only use the hit number ``chosen``
    for the constraints.
    """
    chosen = 0

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