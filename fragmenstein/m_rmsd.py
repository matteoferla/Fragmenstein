from __future__ import annotations
########################################################################################################################

__doc__ = \
    """
Combined RMSD
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################


from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem

from typing import Sequence, List, Optional, Tuple

import json, re

from .monster import Monster

class mRSMD:
    """RMSD are unaligned and in Å.

    The RMSD has been calculated differently.
    The inbuilt RMSD calculations in RDKit (``Chem.rdMolAlign.GetBestRMS``) align the two molecules,
    this does not align them.
    This deals with the case of multiple hits.
    As a comparision, For euclidean distance the square root of the sum of the differences in each coordinates is taken.
    As a comparision, For a regular RMSD the still-squared distance is averaged before taking the root.
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
        currected output of monster.origin_from_mol() or cls.get_origins(to-be-scored-mol, annotated)

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
            tatoms += len(mapping)
            if len(mapping):
                self.rmsds.append((md/ len(mapping)) ** 0.5)
            else:
                self.rmsds.append(None)
        if tatoms:
            self.mrmsd = (sum(mds) / tatoms) ** 0.5
        else:
            self.mrmsd = None


    @classmethod
    def from_unannotated_mols(cls,
                            moved_followup: Chem.Mol,
                            hits: Sequence[Chem.Mol],
                            placed_followup: Chem.Mol
                            ) -> mRSMD:
        """
        Mapping is done by positional overlap between placed_followup and hits
        This mapping is the applied to moved_followup.

        :param moved_followup: The mol to be scored
        :param hits: the hits to score against
        :param placed_followup: the mol to determine how to score
        :return:
        """
        mappings = []
        moved_followup = AllChem.DeleteSubstructs(moved_followup, Chem.MolFromSmiles('*'))
        placed_followup = AllChem.DeleteSubstructs(placed_followup, Chem.MolFromSmiles('*'))
        if moved_followup.GetNumAtoms() != placed_followup.GetNumAtoms():
            # they may differ just because protons
            placed_followup = Chem.AddHs(placed_followup)
            assert moved_followup.GetNumAtoms() == placed_followup.GetNumAtoms(), 'moved and placed are different!'
        for h, hit in enumerate(hits):
            mappings.append(list(Monster.get_positional_mapping(hit, placed_followup).items()))
        return cls(moved_followup, hits, mappings)

    @classmethod
    def from_annotated_mols(cls,
                  annotated_followup: Chem.Mol,
                  hits: Optional[Sequence[Chem.Mol]]=None
                  ) -> mRSMD:
        """
        Monster leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
        This classmethod accepts a followup with has this.

        :param annotated_followup:
        :param hits:
        :return:
        """
        if cls.is_xyz_annotated(annotated_followup):
            return cls.from_internal_xyz(annotated_followup)
        mappings = cls._mapping_from_annotated_and_hits(annotated_followup, hits)
        return cls(annotated_followup, hits, mappings)

    @classmethod
    def is_origin_annotated(cls, mol: Chem.Mol) -> bool:
        for atom in mol.GetAtoms():
            if len(cls._get_origin(atom)) > 0:
                return True
        else:
            return False

    @classmethod
    def is_xyz_annotated(cls, mol: Chem.Mol) -> bool:
        for atom in mol.GetAtoms():
            if len(cls._get_xyz(atom)) > 0:
                return True
        else:
            return False

    @classmethod
    def _mapping_from_annotated_and_hits(cls,
                                         annotated_followup: Chem.Mol,
                                         hits: Sequence[Chem.Mol]):
        assert cls.is_origin_annotated(annotated_followup), 'This molecules is not annotated.'
        mappings = []
        for h, hit in enumerate(hits):
            hname = hit.GetProp('_Name')
            mapping = []
            if hname == '':
                print(f'{hit} has no name!')
            else:
                for i in range(annotated_followup.GetNumAtoms()):
                    atom = annotated_followup.GetAtomWithIdx(i)
                    for oel in cls._get_origin(atom):
                        rex = re.match(hname + '\.(\d+)', oel)
                        if rex is not None:
                            h = int(rex.group(1))
                            mapping.append((i, h))
            mappings.append(mapping)
        return mappings

    @classmethod
    def from_other_annotated_mols(cls,
                            followup: Chem.Mol,
                            hits: Sequence[Chem.Mol],
                            annotated: Chem.Mol
                            ) -> mRSMD:
        # the former way has issues with isomorphisms
        # cls.copy_origins(annotated, followup)
        # return cls.from_annotated_mols(followup, hits)
        targets, _ = cls.copy_all_possible_origins(annotated, followup)
        assert targets, 'Molecule could not be mapped.'
        results = [(target, cls.from_annotated_mols(target, hits)) for target in targets]
        return list(sorted(results, key=lambda x: x[1].mrmsd))[0][1]



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
    def mock(cls):
        self = cls.__new__(cls)
        self.followup = None
        self.hits = None
        self.mappings = {}
        self.mrmsd = float('nan')
        self.rmsds = []
        return self


    @classmethod
    def copy_origins(cls, annotated: Chem.Mol, target: Chem.Mol):
        """
        Monster leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
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
    def copy_all_possible_origins(cls, annotated: Chem.Mol, target: Chem.Mol) -> Tuple[List[Chem.Mol], List[List[int]]]:
        """
        Monster leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
        However, the atom order seems to be maintained but I dont trust it. Also dummy atoms are stripped.

        :param annotated:
        :param target:
        :return: a list of mols and a list of orgins (a list too)
        """
        mcs = rdFMCS.FindMCS([target, annotated],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareAny,
                             ringMatchesRingOnly=True)
        common = Chem.MolFromSmarts(mcs.smartsString)
        options = []
        originss = []
        for target_match in target.GetSubstructMatches(common):
            for anno_match in annotated.GetSubstructMatches(common):
                dmapping = dict(zip(target_match, anno_match))
                origins = []
                option = Chem.Mol(target)
                for i in range(option.GetNumAtoms()):
                    if i in dmapping:
                        atom = annotated.GetAtomWithIdx(dmapping[i])
                        tatom = option.GetAtomWithIdx(i)
                        o = cls._get_origin(atom)
                        tatom.SetProp('_Origin', json.dumps(o))
                        xyz = cls._get_xyz(atom)
                        if xyz:
                            cls._set_xyz(tatom, xyz)
                options.append(option)
                originss.append(origins)
        return options, originss


    @classmethod
    def migrate_origin(cls, mol: Chem.Mol, tag='_Origin') -> Chem.Mol:
        """
        The origin list may be saved as a molecule property rather than an atom -saved as a mol say.

        :param mol: mol to fix
        :param tag: name of prop
        :return: the same mol
        """
        assert mol.HasProp(tag), f'There is no tag {tag}'
        origins = json.loads(mol.GetProp(tag))
        assert len(origins) == mol.GetNumAtoms(), f'Mismatch {len(origins)} vs. {mol.GetNumAtoms()}'
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetProp('_Origin', json.dumps(origins[i]))
        return mol

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

    @classmethod
    def _get_xyz(cls, atom: Chem.Atom) -> Tuple[float]:
        if atom.HasProp('_x'):
            return (atom.GetDoubleProp('_x'),
                    atom.GetDoubleProp('_y'),
                    atom.GetDoubleProp('_z'))
        else:
            return ()

    @classmethod
    def _set_xyz(cls, atom: Chem.Atom, xyz):
        if len(xyz):
            atom.SetDoubleProp('_x', xyz[0]),
            atom.SetDoubleProp('_y', xyz[1]),
            atom.SetDoubleProp('_z', xyz[2])

    @classmethod
    def from_internal_xyz(cls, annotated_followup):
        """
        This is an alternative for when the atoms have _x, _y, _z
        
        :param annotated_followup:
        :return:
        """
        self = cls.__new__(cls)
        self.followup = annotated_followup
        self.hits = []
        self.mappings = []
        self.rmsds = []
        self.mrmsd = float('nan')
        conf = annotated_followup.GetConformer()
        n = 0
        tatoms = 0
        for a, atom in enumerate(annotated_followup.GetAtoms()):
            if atom.HasProp('_x'):
                x, y, z = cls._get_xyz(atom)
                tatoms += 1
                n += sum([(conf.GetAtomPosition(a).x - x) ** 2 +
                            (conf.GetAtomPosition(a).y - y) ** 2 +
                            (conf.GetAtomPosition(a).z - z) ** 2])
        self.mrmsd = (n / tatoms) ** 0.5
        return self


