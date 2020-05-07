from rdkit import Chem
from rdkit.Chem import rdFMCS

from typing import Sequence, List, Optional, Tuple

import json, re

from .core import Fragmenstein

class mRSMD:
    """RMSD are unaligned and in Å.

    The RMSD has been calculated differently.
    The inbuilt RMSD calculations in RDKit (``Chem.rdMolAlign.GetBestRMS``) align the two molecules,
    this does not align them.
    This deals with the case of multiple hits.
    For euclidean distance the square root of the sum of the differences in each coordinates is taken.
    For a regular RMSD the still-squared distance is averaged before taking the root.
    Here the average is done across all the atom pairs between each hit and the followup.
    Therefore, atoms in followup that derive in the blended molecule by multiple atom are scored multiple times.

    $\sqrt{\frac{\sum_{i}^{N_{\rm{hits}}} (\sum_{i}^{n} (q_{i,\rm{x}} - h_{i,\rm{x}})^2 + (q_{i,\rm{y}} - h_{i,\rm{y}})^2 + (q_{i,\rm{z}} - h_{i,\rm{z}})^2 }{n\cdot m}}$



    """

    def __init__(self,
                 followup: Chem.Mol,
                 hits: Sequence[Chem.Mol],
                 mappings: List[List[Tuple[int, int]]]):
        """
        This is not meant to be called directly.
        mappings is a list of len(hits) containing lists of tuples of atom idx that go from followup to hit


        The hit _Name must match that in origin!
        currected output of fragmenstein.origin_from_mol() or cls.get_origins(to-be-scored-mol, annotated)

        :param followup: the followup compounds
        :param hits: the fragment hits
        :param mappings: a complicated affair...
        """
        self.followup = followup
        self.hits = hits
        self.mappings = mappings
        # calculate
        mds = []
        self.rmsds = []
        tatoms = 0
        for hit, mapping in zip(hits, mappings):
            md = self.calculate_msd(self.followup, hit, mapping)
            mds.append(md)
            tatoms = len(mapping)
            self.rmsds.append((md/ len(mapping)) ** 0.5)
        self.mrmsd = (sum(mds)/ tatoms) ** 0.5


    @classmethod
    def from_unnotated_mols(cls,
                            moved_followup: Chem.Mol,
                            hits: Sequence[Chem.Mol],
                            placed_followup: Chem.Mol
                            ):
        mappings = []
        assert moved_followup.GetNumAtoms() == placed_followup.GetNumAtoms(), 'moved and placed are different!'
        for h, hit in enumerate(hits):
            mappings.append(list(Fragmenstein.get_positional_mapping(hit, placed_followup).items()))
        return cls(moved_followup, hits, mappings)

    @classmethod
    def from_annotated_mols(cls,
                  annotated_followup: Chem.Mol,
                  hits: Sequence[Chem.Mol]
                  ):
        """
        Fragmenstein leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
        This classmethod accepts a followup with has this.

        :param annotated_followup:
        :param hits:
        :return:
        """
        mappings = []
        for h, hit in enumerate(hits):
            hname = hit.GetProp('_Name')
            mapping = []
            for i in range(annotated_followup.GetNumAtoms()):
                atom = annotated_followup.GetAtomWithIdx(i)
                for oel in cls._get_origin(atom):
                    if hname in oel:
                        h = int(re.match(hname+'\.(\d+)', oel).group(1))
                        mapping.append((i, h))
            mappings.append(mapping)
        return cls(annotated_followup, hits, mappings)

    @classmethod
    def from_other_annotated_mols(cls,
                            followup: Chem.Mol,
                            hits: Sequence[Chem.Mol],
                            annotated: Chem.Mol
                            ):
        cls.copy_origins(annotated, followup)
        return cls.from_annotated_mols(followup, hits)

    def calculate_msd(self, molA, molB, mapping):
        """
        A nonroot rmsd.

        :param molA:
        :param molB:
        :param mapping: lists of tuples of atom idx that go from molA to molB
        :return: nonroot rmsd
        """
        confA = molA.GetConformer()
        confB = molB.GetConformer()
        return sum([(confA.GetAtomPosition(a).x - confB.GetAtomPosition(b).x) ** 2 +
                    (confA.GetAtomPosition(a).y - confB.GetAtomPosition(b).y) ** 2 +
                    (confA.GetAtomPosition(a).z - confB.GetAtomPosition(b).z) ** 2 for a, b in mapping])

    def calculate_rmsd(self, molA, molB, mapping):
        return (self.calculate_msd(molA, molB, mapping) / len(mapping)) ** 0.5

    @classmethod
    def copy_origins(cls, annotated: Chem.Mol, target: Chem.Mol):
        """
        Fragmenstein leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
        However, the atom order seems to be maintained but I dont trust it. Also dummy atoms are stripped.

        :param annotated:
        :param target:
        :return: a list of origins
        """
        mcs = rdFMCS.FindMCS([target, annotated],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareAny,
                             ringMatchesRingOnly=True)
        common = Chem.MolFromSmarts(mcs.smartsString)
        dmapping = dict(zip(target.GetSubstructMatch(common), annotated.GetSubstructMatch(common)))
        origins = []
        for i in range(target.GetNumAtoms()):
            if i in dmapping:
                atom = annotated.GetAtomWithIdx(dmapping[i])
                tatom = target.GetAtomWithIdx(i)
                o = cls._get_origin(atom)
                tatom.SetProp('_Origin', json.dumps(o))
        return origins

    @classmethod
    def _get_origin(cls, atom: Chem.Atom) -> List[str]:
        if atom.HasProp('_Origin'):
            o = atom.GetProp('_Origin')
            if o != 'none':
                return json.loads(o)
            else:
                return []
        else:
            return []


