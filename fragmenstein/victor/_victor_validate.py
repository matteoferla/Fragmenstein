from rdkit import Chem

from ._victor_base import _VictorBase
from ..m_rmsd import mRSMD


class _VictorValidate(_VictorBase):

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
