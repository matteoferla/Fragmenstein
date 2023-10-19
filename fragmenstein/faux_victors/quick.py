import functools

from rdkit.Chem import rdFMCS
from ..error import FragmensteinError
from ..victor import Victor
from rdkit import Chem, Geometry
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

    @functools.cached_property
    def preminimized_undummied_mol(self) -> Chem.Mol:
        """
        This method is called by the plonking into structure methods.
        Not "positioning" as intended by ``monster`` is done.
        Opening a PDB in RDKit is doable but gets exponentially slow with chain length
        """
        mol = Chem.Mol(self.monster.positioned_mol)
        # ## Store xyz as properties
        atom: Chem.Atom
        conf: Chem.Conformer = mol.GetConformer()  # noqa
        for atom in mol.GetAtoms():  # noqa
            if atom.HasProp('_x'):
                continue
            xyz: Geometry.Point3D = conf.GetAtomPosition( atom.GetIdx() )  # noqa
            atom.SetDoubleProp('_x', float(xyz.x))
            atom.SetDoubleProp('_y', float(xyz.y))
            atom.SetDoubleProp('_z', float(xyz.z))
        # ## Minimise
        if self.monster_mmff_minisation:
            self.journal.debug(f'{self.long_name} - pre-minimising monster (MMFF)')
            min_result = self.monster.mmff_minimize(mol, allow_lax=False, ff_max_displacement=float('nan'))
            if not min_result.success:
                raise FragmensteinError('Could not preminize')
            ddG = min_result.delta
            if ddG > 100:
                raise FragmensteinError('Poor preminization')
        return AllChem.DeleteSubstructs(min_result.mol, Chem.MolFromSmiles('*'))
