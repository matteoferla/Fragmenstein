########################################################################################################################

__doc__ = \
    """
Igor energy minises the blended compound using pyrosetta.
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__citation__ = ""
from ..version import __version__

########################################################################################################################

from typing import Dict, List, Optional, Tuple, Union, Sequence

from ._igor_init import _IgorInit
from .pyrosetta_import import pyrosetta
from ._igor_min import _IgorMin
from ._igor_utils import _IgorUtils, Term

# this contains the init and the two classmethods.


class Igor(_IgorInit, _IgorMin, _IgorUtils):
    """
    Regular Igor(..) accepts pyrosetta pose.
    ``Igor.from_pdbblock(..)`` accepts pdb block as str,
    while ``Igor.from_pdbfile(..)`` accepts filename as str.


    ``ligand`` can be one of many things. default is 'LIG'. But it can be

    * pose index (123)
    * PDB index '123A'
    * a tuple of (PDB resi, PDB chain)
    * a residue name in uppercase "LIG"
    * a pyrosetta.Vector1 where 1 == the ligand.

    If key_residues is None, only the connecting residue is added (if present in the LINK record).
    This is overridden if one of many options are given.
    If it is a pyrosetta.Vector1 it is assumed that 1 mean select this residue (as a result of a ``selector.apply(pose)`` operation)
    If it is a list or tuple, the elements are interpreted similarly to ligand.

    This class inherits the base classes ``_IgorInit``, ``_IgorMin``, ``_IgorUtils``, which contain the initialisation  methods,
    minimisation  methods and extra methods respectively.
    """

    def residues_in_selector(self, pose: pyrosetta.Pose, selector) -> List[str]:
        """
        This method is just for checking purposes for humans basically.
        """
        return self.residues_in_vector(pose, selector.apply(pose))

    def residues_in_vector(self, pose: pyrosetta.Pose, vector: pyrosetta.Vector1) -> List[str]:
        """
        This method is just for checking purposes for humans basically.
        """
        ResidueVector = pyrosetta.rosetta.core.select.residue_selector.ResidueVector
        return [pose.pdb_info().pose2pdb(r) for r in ResidueVector(vector)]
