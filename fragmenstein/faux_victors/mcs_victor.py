from ..victor import Victor
from .mcs_monster import MCSMerger
from rdkit import Chem

class Mictor(Victor):
    uses_pyrosetta = False
    """
    MCS hack.
    This is used for benchmarking in the paper. Not for general use.
    """

    def _calculate_combination_chem(self):
        """
        The rdkit part. Monster is used to combine the hits.
        """
        self.journal.info("modded `_calculate_combination_chem` for `MCSMerger` called")
        self._harmonize_warhead_combine()
        # TODO Does combine not need attachment??
        self.monster = MCSMerger(self.hits)
        self.monster.combine(embed_mode=0)
        # note: self.monster_throw_on_discard check is not applicable here hence absence
        self.post_monster_step()  # empty overridable
        self.mol = self.monster.positioned_mol
        self.mol.SetProp('_Name', self.monster.merged_name)
        self.smiles = Chem.MolToSmiles(self.mol)
        # making folder.
        self.make_output_folder()
        self.unmatched = []
