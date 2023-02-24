from ..victor import Victor
from .mcs_monster import MCSMerger
from rdkit import Chem

class Mictor(Victor):
    # MCS hack

    def _calculate_combination_chem(self):
        """
        The rdkit part. Monster is used to combine the hits.
        """
        print("modded `_calculate_combination_chem` for `MCSMerger` called")
        self._harmonize_warhead_combine()
        # TODO Does combine not need attachment??
        self.monster = MCSMerger(self.hits)
        self.monster.combine(embed_mode=0)
        # self.monster_throw_on_discard is not applicable
        self.mol = self.monster.positioned_mol
        self.mol.SetProp('_Name', self.monster.merged_name)
        self.smiles = Chem.MolToSmiles(self.mol)
        # making folder.
        self.make_output_folder()
        self.unmatched = []
