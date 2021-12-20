########################################################################################################################

__doc__ = \
    """
base methods
    """

########################################################################################################################

import pyrosetta

from typing import Dict, List, Optional, Tuple, Union, Sequence

from warnings import warn


class _IgorInit:
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
        self.pose = pose
        # virtualroot
        if pose.residue(self.pose.total_residue()).name3() != 'XXX':
            pyrosetta.rosetta.core.pose.addVirtualResAsRoot(self.pose)
        pyrosetta.create_score_function('ref2015')(self.pose)
        self.constraint_file = constraint_file
        self.ligand_residue = self._parse_residue(ligand_residue)
        self.key_residues = self._parse_key_residues(key_residues)

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
        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([params_file])
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_file(pose, pdbfile)
        return cls(pose, constraint_file, ligand_residue, key_residues)

    # ============= Private methods for init ===========================================================================

    def _parse_residue(self, residue: Union[int, str, Tuple[int, str], pyrosetta.Vector1]):
        parsed = []
        if residue is None:
            ## assuming it is LIG then.
            ligand_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector()
            ligand_selector.set_residue_name3('LIG')
            m = self._vector2residues(ligand_selector.apply(self.pose))
            if len(m) == 1:
                warn('There are many residues called LIG!. Please specify the name3/resn/resi yourself!')
            for r in m:
                parsed.append(r)
        elif isinstance(residue, int):
            # it is a pose residue index.
            parsed.append(residue)
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
        elif isinstance(residue, pyrosetta.Vector1):
            for r in [i + 1 for i, v in enumerate(residue) if v == 1]:
                parsed.append(r)
        else:
            raise ValueError(f'No idea what {residue} is.')
        parsed = [p for p in parsed if p != 0]
        assert len(parsed), f'There is no {residue} in pose'
        return parsed

    def _parse_key_residues(self,
                            key_residues: Union[None, Sequence[Union[int, str, Tuple[int, str]]],
                                                pyrosetta.Vector1]):
        parsed = []
        if key_residues is None:
            try:
                parsed = [self.pose.residue(self.ligand_residue[0]).connect_map(1).resid()]
            except RuntimeError:
                warn('No covalent bond with the ligand.')
        else:
            for k in key_residues:
                parsed.extend(self._parse_residue(k))
        return parsed
