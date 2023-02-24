from rdkit import Chem
from ..victor import Victor


class Wictor(Victor):
    """
    This Victor does not call Igor

    ... code-block:: python
        from framgenstein.faux_victors import Wictor
        from fragmenstein.demo import Mac1

        Wictor.capture_rdkit_log()
        # Wictor.enable_stdout()
        Wictor.error_to_catch = ()
        wicky = Wictor(hits=Mac1.get_n_filtered_mols(2),
                       pdb_block=Mac1.get_template(),
                       )
        wicky.combine()
    """

    def _calculate_combination_thermo(self):
        # override igor.
        pass

    def _calculate_placement_thermo(self):
        # override igor.
        pass

    def post_monster_step(self):
        # this is a black overridable methods that will be the last thing called
        self.minimized_mol: Chem.Mol = self.monster.positioned_mol
        self.minimized_pdbblock: str = self.unminimized_pdbblock
