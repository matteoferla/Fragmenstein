from rdkit.Chem import rdFMCS
from ..error import FragmensteinError
from ..victor import Victor
from rdkit import Chem
from rdkit.Chem import AllChem


class Quicktor(Victor):
    """
    This class is intended to speed up the Victor placement step, by:

    * doing a strict MCS only
    * not allowing constrained atoms to move during RDKit minimisation and killing them if they fail.
    * doing a quick pyrosetta minimisation

    Example usage:

    .. code-block:: python
       lab = Laboratory(hits, pdb_block=pdb_block)
       lab.Victor = Quicktor
       to_be_placed: pd.DataFrame = place_input_validator(to_be_placed)
       placed: pd.DataFrame = lab.place(to_be_placed_df, n_cores=55, expand_isomers=False, timeout=60)
    """
    quick_reanimation = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.monster.matching_modes = [{'atomCompare': rdFMCS.AtomCompare.CompareElements,
                                        'bondCompare': rdFMCS.BondCompare.CompareOrder,
                                        'ringCompare': rdFMCS.RingCompare.PermissiveRingFusion,
                                        'ringMatchesRingOnly': True}]
        self.monster_mmff_minisation = True  # as is default. Repeated for emphasis.


    def _get_preminimized_undummied_monster(self) -> Chem.Mol:
        """
        This method is called by the plonking into structure methods.
        Not "positioning" as intended by ``monster`` is done.
        Opening a PDB in RDKit is doable but gets exponentially slow with chain length
        """
        mol = Chem.Mol(self.monster.positioned_mol)
        if self.monster_mmff_minisation:
            self.journal.debug(f'{self.long_name} - pre-minimising monster (MMFF)')
            if not self.monster.mmff_minimize(mol, allow_lax=False, ff_dist_thr=float('nan')):
                raise FragmensteinError('Could not preminize')
            ddG = self.monster.MMFF_score(mol, delta=True)
            if ddG > 100:
                raise FragmensteinError('Poor preminization')
        return AllChem.DeleteSubstructs(mol, Chem.MolFromSmiles('*'))
