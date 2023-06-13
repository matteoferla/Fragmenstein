########################################################################################################################

__doc__ = \
    """
base methods
    """

########################################################################################################################
import os.path

from .pyrosetta_import import pyrosetta  # the real mcCoy or a mock.
from ._igor_base import _IgorBase

from typing import Dict, List, Optional, Tuple, Union, Sequence

from warnings import warn


class _IgorInit(_IgorBase):
    atom_pair_constraint = 10
    angle_constraint = 10
    coordinate_constraint = 1
    fa_intra_rep = 0.005

    # ============= Init ===============================================================================================

    def __init__(self,
                 pose: pyrosetta.Pose,
                 constraint_file: str,
                 ligand_residue: Union[str, int, Tuple[int, str], pyrosetta.Vector1] = 'LIG',
                 key_residues: Union[None, Sequence[Union[int, str, Tuple[int, str]]], pyrosetta.Vector1] = None):
        """
        Given a pose with a blended ligand at ligand residue. Load it (ready for minimisation).

        :param pose: pose.
        :param constraint_file: filename
        :param ligand_residue: ligand -see class docstring
        :param key_residues: multiple entries -see class docstring
        """
        self.pose = pose  #: pyrosetta.Pose
        # virtualroot
        if pose.residue(self.pose.total_residue()).name3() != 'XXX':
            pyrosetta.rosetta.core.pose.addVirtualResAsRoot(self.pose)
        pyrosetta.create_score_function('ref2015')(self.pose)
        self.constraint_file = constraint_file  #: str
        self.ligand_residue: List[int] = self._parse_residue(ligand_residue)
        self.key_residues: List[int] = self._parse_key_residues(key_residues)

    @classmethod
    def from_pdbblock(cls,
                      pdbblock: str,
                      params_file: str,
                      constraint_file: str,
                      ligand_residue: Union[str, int, Tuple[int, str], pyrosetta.Vector1] = 'LIG',
                      key_residues: Union[None, Sequence[Union[int, str, Tuple[int, str]]], pyrosetta.Vector1] = None):


        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([params_file])
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, pdbblock)
        return cls(pose, constraint_file, ligand_residue, key_residues)

    @classmethod
    def from_pdbfile(cls,
                     pdbfile: str,
                     params_file: str,
                     constraint_file: str,
                     ligand_residue: Union[str, int, Tuple[int, str], pyrosetta.Vector1] = 'LIG',
                     key_residues: Union[None, Sequence[Union[int, str, Tuple[int, str]]], pyrosetta.Vector1] = None):
        """
        Given a PDB with the ligand load it.

        :param pdbfile: pdb file
        :param params_file: params file
        :param constraint_file: filename
        :param ligand_residue: ligand -see class docstring
        :param key_residues: multiple entries -see class docstring
        :return:
        """
        for filename in (pdbfile, constraint_file, params_file):
            if not os.path.exists(filename):
                raise FileNotFoundError(f'{filename} does not exist')
            if os.stat(filename).st_size < 10:
                raise ValueError(f'{filename} is only {os.stat(filename).st_size} bytes long')
        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([params_file])
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_file(pose, pdbfile)
        return cls(pose, constraint_file, ligand_residue, key_residues)

    # ============= Private methods for init ===========================================================================

    def _parse_residue(self, residue: Union[int, str, Tuple[int, str], pyrosetta.Vector1]):
        parsed: List[int] = []
        if residue is None:
            ## assuming it is LIG then.
            ligand_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector()
            ligand_selector.set_residue_name3('LIG')
            m = self._vector2residues(ligand_selector.apply(self.pose))  # noqa its in Igor_min
            if len(m) > 1:
                warn('There are many residues called LIG!. Please consider specifying the name3/resn/resi yourself!')
            for r in m:
                parsed.append(r)
        elif isinstance(residue, int):
            # it is a pose residue index.
            parsed.append(residue)
        elif isinstance(residue, str) and residue.isdigit():
            parsed.append(int(residue))
        elif isinstance(residue, str) and residue[:-1].isdigit():
            # it is a pdb residue index + chain e.g. 23A.
            res = int(residue[:-1])
            chain = residue[-1]
            r = self.pose.pdb_info().pdb2pose(res=res, chain=chain)
            parsed.append(r)
        elif isinstance(residue, str) and len(residue) in (3, 4) and residue.upper() == residue:
            # it is a residue name, such as LIG
            ligand_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector()
            ligand_selector.set_residue_name3(residue)
            for r in [i + 1 for i, v in enumerate(ligand_selector.apply(self.pose)) if v == 1]:
                parsed.append(r)
        elif isinstance(residue, tuple) and len(residue) == 2:
            parsed = self.pose.pdb_info().pdb2pose(res=residue[0], chain=residue[1])
        # elif isinstance(residue, pyrosetta.Vector1): will raise an error
        # as the latter is a function. that redirects to `vector1_int`, `vector1_double` etc.
        elif hasattr(residue, '__class__') and 'vector1' in residue.__class__.__name__:
            for r in [int(i + 1) for i, v in enumerate(residue) if int(v) == 1]:
                parsed.append(r)
        else:
            raise ValueError(f'No idea what {residue} is.')
        parsed = [p for p in parsed if p != 0]
        assert len(parsed), f'There is no {residue} in pose'
        return parsed

    def _parse_key_residues(self,
                            key_residues: Union[None, Sequence[Union[int, str, Tuple[int, str]]],
                                                pyrosetta.Vector1]):
        parsed: List[int] = []
        residue = self.pose.residue(self.ligand_residue[0])
        if key_residues is None and residue.n_current_residue_connections():
            parsed = [residue.connect_map(1).resid()]
        elif key_residues is None:
            parsed = [1]
        else:
            for k in key_residues:
                parsed.extend(self._parse_residue(k))
        return parsed

