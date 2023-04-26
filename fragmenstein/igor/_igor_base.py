# this is to keep mypy happy
from .pyrosetta_import import pyrosetta
import warnings
from typing import List, Tuple

class _IgorBase:
    def __init__(self):
        warnings.warn('Virtual only method: THIS METHOD SHOULD NOT BE RUN. CALL `Igor`', category=SyntaxWarning)
        self.pose = pyrosetta.Pose()
        self.constraint_file = ''
        self.ligand_residue: List[int] = []
        self.key_residues: List[int] = []
        self.atom_pair_constraint = 10
        self.angle_constraint = 10
        self.coordinate_constraint = 1
        self.fa_intra_rep = 0.005 # default

    def _vector2residues(self, vector: pyrosetta.Vector1) -> List[int]:
        """
        This method is for machines. See ``residues_in_vector`` instead.
        """
        return [i + 1 for i, v in enumerate(vector) if v == 1]