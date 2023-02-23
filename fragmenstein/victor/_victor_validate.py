from rdkit import Chem
from rdkit.Chem import AllChem

from ._victor_base import _VictorBase
from ..m_rmsd import mRMSD
from typing import *


class _VictorValidate(_VictorBase):

    def validate(self, reference_mol: Chem.Mol) -> Dict[str, float]:
        """
        Get how well the results compare.
        Alternative, do a docking with victor.dock() (-> Chem.Mol)

        :param reference_mol: Crystal structure mol
        :return:
        """
        self.reference_mol : Chem.Mol = reference_mol
        return dict(minimized2hits_rmsd=self._get_minimized2hits_rmsd(),
                    reference2hits_rmsd=self._get_reference2hits_rmsd(),
                    reference2minimized_rmsd=self._get_reference2minimized_rmsd(),
                    )

    def _calc_rmsd_safeguard(self, func: Callable, *args, **kwargs):
        try:
            # compare with reference mol
            return func(*args, **kwargs).mrmsd
        except KeyboardInterrupt as err:
            raise err
        except self.error_to_catch as err:
            self.journal.error(f'{err.__class__.__name__}: {err} in validation step.')
            pass
        return float('nan')

    def _get_reference2hits_rmsd(self) -> float:
        return self._calc_rmsd_safeguard(mRMSD.from_other_annotated_mols,
                                         followup=self.reference_mol,
                                         hits=self.hits,
                                         annotated=self.monster.positioned_mol)

    def _get_reference2minimized_rmsd(self) -> float:
        return self._calc_rmsd_safeguard(mRMSD.from_unannotated_same_mols,
                                         moved_followup=self.minimized_mol,
                                         placed_followup=self.reference_mol)

    def _get_minimized2hits_rmsd(self) -> float:
        # actually is annotated, but best do the full loop in case it was tinkered with.
        return self._calc_rmsd_safeguard(mRMSD.from_other_annotated_mols,
                                         followup=self.minimized_mol,
                                         hits=self.hits,
                                         annotated=self.monster.positioned_mol)


