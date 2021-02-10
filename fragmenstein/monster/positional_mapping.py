########################################################################################################################

__doc__ = \
    """
Positional mapping
    """

########################################################################################################################


from rdkit import Chem
from typing import Dict, List, Tuple
import numpy as np


# ========= Get positional mapping =================================================================================

class GPM:
    """
    This class simply contains ``get_positional_mapping`` and is inherited both by Monster and Unmerge.
    ``get_positional_mapping`` teturns a map to convert overlapping atom of A onto B
    """

    cutoff = 2

    @classmethod
    def get_positional_mapping(cls, mol_A: Chem.Mol, mol_B: Chem.Mol, dummy_w_dummy=True) -> Dict[int, int]:
        """
        Returns a map to convert overlapping atom of A onto B
        Cutoff 2 &Aring; (see class attr.)

        :param mol_A: first molecule (Chem.Mol) will form keys
        :param mol_B: second molecule (Chem.Mol) will form values
        :param dummy_w_dummy: match */R with */R.
        :return: dictionary mol A atom idx -> mol B atom idx.
        """

        mols = [mol_A, mol_B]
        confs = [m.GetConformers()[0] for m in mols]
        distance_protomatrix = []
        dummy_distance_protomatrix = []
        ring_distance_protomatrix = []
        for i in range(mols[0].GetNumAtoms()):
            distance_protovector = []
            dummy_distance_protovector = []
            ring_distance_protovector = []
            for j in range(mols[1].GetNumAtoms()):
                distance, dummy_distance, ring_distance = cls._gpm_distance(mols, confs, i, j, dummy_w_dummy)
                distance_protovector.append(distance)
                dummy_distance_protovector.append(dummy_distance)
                ring_distance_protovector.append(ring_distance)
            distance_protomatrix.append(distance_protovector)
            dummy_distance_protomatrix.append(dummy_distance_protovector)
            ring_distance_protomatrix.append(ring_distance_protovector)
        distance_matrix = np.array(distance_protomatrix)
        dummy_distance_matrix = np.array(dummy_distance_protomatrix)
        ring_distance_matrix = np.array(ring_distance_protomatrix)
        if dummy_w_dummy:
            return {**cls._gpm_covert(distance_matrix, cls.cutoff),
                    **cls._gpm_covert(dummy_distance_matrix, cls.cutoff * 2),
                    **cls._gpm_covert(ring_distance_matrix, cls.cutoff)}
        else:
            return {**cls._gpm_covert(distance_matrix, cls.cutoff),
                    **cls._gpm_covert(ring_distance_matrix, cls.cutoff)}

    @classmethod
    def _gpm_distance(cls, mols: List[Chem.Mol], confs: [Chem.Conformer], i, j, dummy_w_dummy=True) \
            -> Tuple[float, float, float]:
        """
        See get_positional_distance

        :param mols:
        :param dummy_w_dummy:
        :return:
        """
        measure_distance = lambda a, b: ((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2) ** 0.5
        is_collapsed_ring = lambda atom: atom.HasProp('_ori_i') and atom.GetIntProp('_ori_i') == -1
        zeroth = mols[0].GetAtomWithIdx(i)
        first = mols[1].GetAtomWithIdx(j)
        if dummy_w_dummy and \
                zeroth.GetSymbol() == '*' and \
                first.GetSymbol() == '*':
            distance = 9999
            dummy_distance = measure_distance(confs[0].GetAtomPosition(i), confs[1].GetAtomPosition(j))
            ring_distance = 9999
        elif dummy_w_dummy and \
                (zeroth.GetSymbol() == '*' or
                 first.GetSymbol() == '*'):
            distance = 9999
            dummy_distance = 9999
            ring_distance = 9999
        elif is_collapsed_ring(zeroth) and is_collapsed_ring(first):
            distance = 9999
            dummy_distance = 9999
            ring_distance = measure_distance(confs[0].GetAtomPosition(i), confs[1].GetAtomPosition(j))
        elif is_collapsed_ring(zeroth) or is_collapsed_ring(first):
            distance = 9999
            dummy_distance = 9999
            ring_distance = 9999
        else:
            distance = measure_distance(confs[0].GetAtomPosition(i), confs[1].GetAtomPosition(j))
            dummy_distance = 9999
            ring_distance = 9999
        return distance, dummy_distance, ring_distance

    @classmethod
    def _gpm_covert(cls, array: np.array, cutoff: float) -> Dict[int, int]:
        """
        See get_positional_distance

        :param array:
        :param cutoff:
        :return:
        """
        # find the closest
        mapping = {}
        while 1 == 1:
            d = np.nanmin(array)
            if d > cutoff:
                break
            w = np.where(array == d)
            f, s = w[0][0], w[1][0]
            mapping[int(f)] = int(s)  # np.int64 --> int
            array[f, :] = np.ones(array.shape[1]) * 999
            array[:, s] = np.ones(array.shape[0]) * 999
        return mapping
