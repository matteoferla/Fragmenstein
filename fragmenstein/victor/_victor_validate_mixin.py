from ._victor_base_mixin import _VictorBaseMixin
from ..monster import Monster
from ..igor import Igor
from ..m_rmsd import mRSMD
from ..rectifier import Rectifier
from typing import List, Optional, Dict, Union, Callable
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_to_params import Params, Constraints
import time, warnings, os

class _VictorValidateMixin(_VictorBaseMixin):

    def validate(self, reference_mol: Chem.Mol):
        """
        Get how well the results compare.
        Alternative, do a docking with victor.dock() (-> Chem.Mol)

        :param reference_mol: Crystal structure mol
        :return:
        """
        try:
            # compare with reference mol
            return mRSMD.from_other_annotated_mols(reference_mol,
                                                  self.hits,
                                                  self.monster.positioned_mol).mrmsd
        except self.error_to_catch as err:
            self.journal.error(f'{err.__class__.__name__}: {err} in validation step.')
            pass
        return float('nan')
