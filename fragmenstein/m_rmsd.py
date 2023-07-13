from __future__ import annotations
from collections import defaultdict
########################################################################################################################

__doc__ = \
    """
Combined RMSD
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__citation__ = ""

import warnings

from .version import __version__

########################################################################################################################


from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem
from rdkit.Geometry import Point3D

from typing import Sequence, List, Optional, Tuple

import json, re
import numpy as np

from .monster import Monster
from typing import *

class mRMSD:
    """Mols are non-superposed ('aligned') for the RMSD and in Å.

    The RMSD has been calculated differently.
    The inbuilt RMSD calculations in RDKit (``Chem.rdMolAlign.GetBestRMS``) align the two molecules,
    this does not align/superpose them.
    This deals with the case of multiple hits.
    As a comparision, For euclidean distance the square root of the sum of the differences in each coordinates is taken.
    As a comparision, For a regular RMSD the still-squared distance is averaged before taking the root.
    Here the average is done across all the atom pairs between each hit and the followup.
    Therefore, atoms in followup that derive in the blended molecule by multiple atom are scored multiple times.

    .. math::

        \sqrt{\frac{\sum_{i}^{N_{\rm{hits}}} (\sum_{i}^{n} (q_{i,\rm{x}} - h_{i,\rm{x}})^2 +  \\
        (q_{i,\rm{y}} - h_{i,\rm{y}})^2 + (q_{i,\rm{z}} - h_{i,\rm{z}})^2 }{n\cdot m}}



    """
    # faux enum:
    HIT_BASED = 0  # regularly instantiated
    XYZ_BASED = 1  # from_internal_xyz
    IDENTITY_BASED = 2  # identity

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
        self.mode = self.HIT_BASED
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

    @staticmethod
    def generate_overlap_mapping(mol_a: Chem.Mol, mol_b: Chem.Mol) -> List[Tuple[int, int]]:
        # Monster.get_positional_mapping : Dict[int, int]
        map_dict : Dict[int, int] = Monster.get_positional_mapping(mol_a, mol_b)
        return list(map_dict.items())

    @classmethod
    def from_unannotated_mols(cls,
                            moved_followup: Chem.Mol,
                            hits: Sequence[Chem.Mol],
                            placed_followup: Chem.Mol
                            ) -> mRMSD:
        """
        Mapping is done by positional overlap between placed_followup and hits
        This mapping is the applied to moved_followup.
        The mapping is not between placed and moved. But the former acts as a go between.

        :param moved_followup: The mol to be scored
        :param hits: the hits to score against
        :param placed_followup: the mol to determine how to score
        :return:
        """
        placed2moved : Dict[int, int] = dict(cls.from_unannotated_same_mols(placed_followup=placed_followup,
                                                                            moved_followup=moved_followup
                                                                            ).mappings[0]
                                             )
        mappings = []
        for h, hit in enumerate(hits):
            placed2hit : List[Tuple[int, int]] = cls.generate_overlap_mapping(placed_followup, hit)
            mappings.append([(placed2moved[pi], hi) for pi, hi in placed2hit if pi in placed2moved])
        return cls(followup=moved_followup, hits=hits, mappings=mappings)

    @classmethod
    def from_unannotated_same_mols(cls,
                            moved_followup: Chem.Mol,
                            placed_followup: Chem.Mol
                            ) -> mRMSD:
        """
        Mapping is not done by positional overlap between placed_followup and moved_followup
        But by shape overlap that yields the lowest RMSD
        This means that every atom matches reguardless of hits.

        :param moved_followup: The mol to be scored
        :param placed_followup: the mol to determine how to score
        :return:
        """
        # placed idx -> moved idx
        munge_mapping = lambda match: list(zip(range(placed_followup.GetNumAtoms()), match))
        candidate_mappings = [munge_mapping(match) for match in moved_followup.GetSubstructMatches(placed_followup)]
        assert candidate_mappings, 'No matches?!'
        # mappings is places to moved hence the order
        candidate_mrmsds = [cls(followup=placed_followup, hits=[moved_followup], mappings=[mapping])
                            for mapping in candidate_mappings]
        self = sorted(candidate_mrmsds, key=lambda m: m.mrmsd)[0]
        # self.hits is placed_followup not hits
        self.mode = self.IDENTITY_BASED
        return self

    @classmethod
    def from_annotated_mols(cls,
                  annotated_followup: Chem.Mol,
                  hits: Optional[Sequence[Chem.Mol]]=None
                  ) -> mRMSD:
        """
        Monster leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
        This classmethod accepts a followup with has this.

        :param annotated_followup:
        :param hits:
        :return:
        """
        if cls.is_xyz_annotated(annotated_followup):
            self = cls.from_internal_xyz(annotated_followup)
            self.hits = hits
            return self
        mappings = cls._mapping_from_annotated_and_hits(annotated_followup, hits)
        return cls(annotated_followup, hits, mappings)

    @classmethod
    def is_origin_annotated(cls, mol: Chem.Mol) -> bool:
        """
        This is atom not mol.
        """
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
                            annotated: Chem.Mol,
                            ) -> mRMSD:
        """
        The two molecules are the same (atom names included) but the followup lacks annotations.
        """
        cls.overannotate(annotated, hits=hits, priority=('atom_origin', 'mol_origin', 'xyz'))
        if all([cls.is_xyz_annotated(annotated),
                annotated.GetAtomWithIdx(0).GetPDBResidueInfo() is not None,
                followup.GetAtomWithIdx(0).GetPDBResidueInfo() is not None
                ]):
            # match by atom name
            name2xyz = {}
            for atom in annotated.GetAtoms():
                pi = atom.GetPDBResidueInfo()
                if pi and atom.HasProp('_x'):
                    name2xyz[pi.GetName()] = {cart: atom.GetDoubleProp("_" + cart) for cart in ('x', 'y', 'z')}
            for atom in followup.GetAtoms():
                pi = atom.GetPDBResidueInfo()
                if pi is None:
                    continue
                name = pi.GetName()
                if pi and name in name2xyz:
                    for cart in ('x', 'y', 'z'):
                        atom.SetDoubleProp('_' + cart, name2xyz[name][cart])
            return cls.from_annotated_mols(followup, hits)
        elif (atom.GetSymbol() for atom in annotated.GetAtoms()) == (atom.GetSymbol() for atom in followup.GetAtoms()):
            # element sequence matches... but has no atom names
            for anno_atom, follow_atom in zip(annotated.GetAtoms(), followup.GetAtoms()):
                if anno_atom.HasProp('_Origin'):
                    follow_atom.SetProp('_Origin', anno_atom.GetProp('_Origin'))
                if anno_atom.HasProp('_x'):
                    follow_atom.SetProp('_x', anno_atom.GetProp('_x'))
                    follow_atom.SetProp('_y', anno_atom.GetProp('_y'))
                    follow_atom.SetProp('_z', anno_atom.GetProp('_z'))
        else:
            # has become nearly redundant with from_unannotated_mols except this expect annotated to have:
            # from_unannotated_mols
            # the former way has issues with isomorphisms
            # cls.copy_origins(annotated, followup)
            # return cls.from_annotated_mols(followup, hits)
            posibilities : Tuple[List[Chem.Mol], List[List[int]]]= cls.copy_all_possible_origins(annotated, followup)
            targets: List[Chem.Mol] = posibilities[0]
            assert targets, 'Molecule could not be mapped.'
            results = [(target, cls.from_annotated_mols(target, hits)) for target in targets]
            return list(sorted(results, key=lambda x: x[1].mrmsd))[0][1]


    def calculate_msd(self, molA, molB, mapping) -> float:
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

    def calculate_rmsd(self, molA, molB, mapping) -> float:
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
        > This is no longer used by ``from_other_annotated_mols``.

        Monster leaves a note of what it did. atom prop _Origin is a json of a list of mol _Name dot AtomIdx.
        However, the atom order seems to be maintained but I dont trust it. Also dummy atoms are stripped.

        :param annotated:
        :param target:
        :return: a list of origins
        """
        warnings.warn('This is no longer used by ``from_other_annotated_mols``.', DeprecationWarning)
        mcs = rdFMCS.FindMCS([target, annotated],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareAny,
                             ringMatchesRingOnly=True)
        common = Chem.MolFromSmarts(mcs.smartsString)
        dmapping = dict(zip(target.GetSubstructMatch(common), annotated.GetSubstructMatch(common)))
        origins = []
        cls.overannotate(annotated, [], ('atom_origin', 'mol_origin', 'surrender'))
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
        cls.overannotate(annotated, [], ('atom_origin', 'mol_origin', 'xyz'))
        mcs = rdFMCS.FindMCS([target, annotated],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareAny,
                             ringMatchesRingOnly=True)
        common = Chem.MolFromSmarts(mcs.smartsString)
        options = []
        originss = []
        # adding Unique=False will make benzene isomorphism not an issue, but the code will grind to a halt.
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
    def overannotate(cls, mol: Chem.Mol, hits: List[Chem.Mol], priority=('mol_origin', 'atom_origin', 'xyz')):
        """
        Unfortunately, in an attempt to make users happy, I have to make the code more complicated.
        There are three different annotation systems going on.
        The first is the property '_Origin' in the ``Chem.Mol``.
        The second is the property '_Origin' in the ``Chem.Atom``.
        The third is the properties '_x', '_y', '_z' in the ``Chem.Atom``.

        This method will take the first priority that is available and propagate it to the other two.

        :param mol: the annotated molecule
        :param hits: the hits that are used to annotate the molecule
        :param priority: the priority of the annotation, a sequence of up to four strings:
            'mol_origin', 'atom_origin', 'xyz', 'surrender'
        """
        if len(priority) == 0:
            raise ValueError('molecule is not annotated')
        task = priority[0]
        # -----------------------------------------------
        if task == 'mol_origin':
            if not mol.HasProp('_Origin'):
                return cls.overannotate(mol, hits, tuple(priority[1:]))
            origins = json.loads(mol.GetProp('_Origin'))
            assert len(origins) == mol.GetNumAtoms(), f'Mismatch {len(origins)} vs. {mol.GetNumAtoms()}'
            for i, atom in enumerate(mol.GetAtoms()):
                # add origin to atom
                atom.SetProp('_Origin', json.dumps(origins[i]))
                # add xyz to atom
                cls._atom_origin_to_xyz(origins[i], atom, hits)
            return None
        # -----------------------------------------------
        elif task == 'atom_origin':
            if not cls.is_origin_annotated(mol):
                return cls.overannotate(mol, hits, tuple(priority[1:]))
            origins = [cls._get_origin(atom) for atom in mol.GetAtoms()]
            # add origin to mol
            mol.SetProp('_Origin', json.dumps(origins))
            # add xyz to atom
            for atom, origin in zip(mol.GetAtoms(), origins):
                cls._atom_origin_to_xyz(origin, atom, hits)
        # -----------------------------------------------
        elif task == 'xyz':
            if not cls.is_xyz_annotated(mol):
                return cls.overannotate(mol, hits, tuple(priority[1:]))
            # this cannot be consolidated...
            pass
        elif task == 'surrender':
            return None
        else:
            raise ValueError(f'Unknown task {task}')

    @classmethod
    def _atom_origin_to_xyz(cls, origins: List[str], atom: Chem.Atom, hits: List[Chem.Mol]):
        points = []
        if len(hits) == 0:
            return  # no hits
        hitdex: Dict[str, Chem.Mol] = {hit.GetProp('_Name'): hit for hit in hits}
        for origin in origins:
            if origin == 'none':
                # thank you human.
                continue
            hit_name, idx = re.match(r'(.*)\.(\d+)', origin).groups()
            atom.SetProp('_ori_name', hit_name)
            if hit_name not in hitdex:
                return
            ref_conf = hitdex[hit_name].GetConformer()
            ref_xyz: Point3D = ref_conf.GetAtomPosition(int(idx))
            points.append(list(ref_xyz))  # noqa this is okay
        if points:
            xyz = np.mean(points, axis=0)
            cls._set_xyz(atom, xyz)

    @classmethod
    def migrate_origin(cls, mol: Chem.Mol, tag='_Origin') -> Chem.Mol:
        """
        The origin list may be saved as a molecule property rather than an atom -saved as a mol say.

        :param mol: mol to fix
        :param tag: name of prop
        :return: the same mol
        """
        raise NotImplementedError('This method was deprecated. Use ``overannotate`` instead. Requires HITS!')

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
    def _get_xyz(cls, atom: Chem.Atom) -> Tuple[float,float,float]:
        if atom.HasProp('_x'):
            return (atom.GetDoubleProp('_x'),
                    atom.GetDoubleProp('_y'),
                    atom.GetDoubleProp('_z'))
        else:
            return ()

    @classmethod
    def _set_xyz(cls, atom: Chem.Atom, xyz):
        if len(xyz):
            atom.SetDoubleProp('_x', float(xyz[0])),
            atom.SetDoubleProp('_y', float(xyz[1])),
            atom.SetDoubleProp('_z', float(xyz[2]))

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
        self.mode = self.XYZ_BASED
        conf = annotated_followup.GetConformer()
        n = 0
        tatoms = 0
        ns = defaultdict(int)
        ts = defaultdict(int)
        for a, atom in enumerate(annotated_followup.GetAtoms()):
            if not atom.HasProp('_x'):
                continue
            x, y, z = cls._get_xyz(atom)
            tatoms += 1
            d = sum([(conf.GetAtomPosition(a).x - x) ** 2 +
                        (conf.GetAtomPosition(a).y - y) ** 2 +
                        (conf.GetAtomPosition(a).z - z) ** 2])
            n += d
            if atom.HasProp('_Origin'):
                origins = json.loads(atom.GetProp('_Origin'))
                for origin in origins:
                    if origin == 'none':
                        continue
                    hit_name, idx = re.match(r'(.*)\.(\d+)', origin).groups()
                    ns[hit_name] += d
                    ts[hit_name] += 1
            elif atom.HasProp('_ori_name'):
                hit_name = atom.GetProp('_ori_name')
                ns[hit_name] += d
                ts[hit_name] += 1
        self.mrmsd = (n / tatoms) ** 0.5
        self.rmsds = [(ns[hit_name] / ts[hit_name]) ** 0.5 for hit_name in ns]
        return self


