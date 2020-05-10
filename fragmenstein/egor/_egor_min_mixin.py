from typing import Dict, List, Optional, Tuple, Union, Sequence

from warnings import warn

import pyrosetta

from rdkit import Chem
from rdkit.Chem import AllChem


class _EgorMinMixin:

    def __init__(self):
        warn('THIS METHOD SHOULD NOT BE RUN. INHERIT _init_')
        self.pose = pyrosetta.Pose()
        self.constraint_file = ''
        self.ligand_residue = []
        self.key_residues = []
        self.atom_pair_constraint = 10
        self.angle_constraint = 10
        self.coordinate_constraint = 1

    def pose2str(self, pose: Optional[pyrosetta.Pose] = None) -> str:
        """
        Convert a pose to a string. Why isn't this standard?

        :param pose: if no pose is provided self.pose is used.
        :return: pdb block
        :rtype: str
        """
        if pose is None:
            pose = self.pose
        buffer = pyrosetta.rosetta.std.stringbuf()
        pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
        return buffer.str()

    def _vector2residues(self, vector: pyrosetta.Vector1) -> List[int]:
        """
        This method is for machines. See ``residues_in_vector`` instead.
        """
        return [i + 1 for i, v in enumerate(vector) if v == 1]

    def mol_from_pose(self) -> Chem.Mol:
        """

        :return: ligand
        :rtype: Chem.Mol
        """
        # mod = self.make_ligand_only_pose()
        # pdbblock = self.pose2str(mod)
        # mol = Chem.MolFromPDBBlock(pdbblock, proximityBonding=False)
        # assert mol is not None, 'Molecule too horrendous to load.'
        # return mol
        mol = Chem.MolFromPDBBlock(self.pose2str(), proximityBonding=False, removeHs=False)
        name3 = self.pose.residue(self.ligand_residue[0]).name3()
        return Chem.SplitMolByPDBResidues(mol, whiteList=[name3])[name3]

    def make_ligand_only_pose(self) -> pyrosetta.Pose:
        """

        :return:
        """
        # this is really janky. But I could not manage to do this without segfaults!
        raise NotImplementedError('This sporadically does not work.')
        mod = self.pose.clone()
        name3 = self.pose.residue(self.ligand_residue[0]).name3()
        ligand_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector(name3)
        not_selector = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(ligand_selector)
        residues = self._vector2residues(not_selector.apply(mod))
        while residues:
            print(residues)
            rbegin = residues[0]
            prev = rbegin
            i = 1
            print(len(residues),i+1 ,prev + 1 , residues[i])
            while len(residues) > i and prev + 1 == residues[i]:
                prev = residues[i]
                i += 1
            rend = residues[i - 1]
            #mod.delete_residue_slow(r)
            mod.delete_residue_range_slow(rbegin, rend)
            residues = self._vector2residues(not_selector.apply(mod))
        return mod

    def MMFF_score(self, mol: Optional[Chem.Mol] = None, delta: bool = False) -> float:
        """
        Merck force field. Chosen over Universal for no reason at all.

        :param mol: ligand
        :type mol: Chem.Mol optional. If absent extracts from pose.
        :param delta: report difference from unbound (minimized)
        :type delta: bool
        :return: kcal/mol
        :rtype: float
        """
        if mol is None:
            mol = self.mol_from_pose()
        AllChem.UFFGetMoleculeForceField(mol)
        ff = AllChem.UFFGetMoleculeForceField(mol)
        ff.Initialize()
        # print(f'MMFF: {ff.CalcEnergy()} kcal/mol')
        if delta:
            pre = ff.CalcEnergy()
            ff.Minimize()
            post = ff.CalcEnergy()
            return pre - post
        else:
            return ff.CalcEnergy()

    def _get_scorefxn(self, name: str = "ref2015"):
        scorefxn = pyrosetta.create_score_function(name)
        # ref2015_cart_cst.wts
        # constrain
        if self.constraint_file:
            self.pose.dump_pdb('test.pdb')
            setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
            setup.constraint_file(self.constraint_file)
            setup.apply(self.pose)
        coord = pyrosetta.rosetta.protocols.relax.AtomCoordinateCstMover()
        coord.set_native_pose(self.pose.clone())
        coord.apply(self.pose)
        stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
        scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), self.atom_pair_constraint)
        scorefxn.set_weight(stm.score_type_from_name("angle_constraint"), self.angle_constraint)
        scorefxn.set_weight(stm.score_type_from_name("coordinate_constraint"), self.coordinate_constraint)
        return scorefxn

    def _get_selector(self,
                      ligand_only: bool = False) -> pyrosetta.rosetta.core.select.residue_selector.ResidueSelector:
        """
        :param ligand_only: selector for ligand only, or with key residues too?
        :return: Gets a selector corresponding to the key_residues
        """
        selectors = []
        if ligand_only:
            it = self.ligand_residue
        else:
            it = self.key_residues + self.ligand_residue
        for pose_residue in it:
            r_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
            r_selector.set_index(pose_residue)
            selectors.append(r_selector)
        or_selector = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector()
        for r_selector in selectors:
            or_selector.add_residue_selector(r_selector)
        return or_selector

    def _get_movemap(self) -> pyrosetta.MoveMap:
        movemap = pyrosetta.MoveMap()
        selector = self._get_selector()
        x = selector.apply(self.pose)
        movemap.set_bb(allow_bb=x)
        movemap.set_chi(allow_chi=x)
        movemap.set_jump(True)
        return movemap

    def _get_movemapfactory(self) -> pyrosetta.rosetta.core.select.movemap.MoveMapFactory:
        selector = self._get_selector()
        mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
        true = pyrosetta.rosetta.core.select.movemap.move_map_action(True)
        mmf.all_chi(False)
        mmf.all_bb(False)
        mmf.all_bondangles(False)
        mmf.all_bondlengths(False)
        mmf.add_chi_action(true, selector)
        mmf.add_bb_action(true, selector)
        mmf.add_bondangles_action(true, selector)
        mmf.add_bondlengths_action(true, selector)
        mmf.set_cartesian(True)
        mmf.add_jump_action(true, pyrosetta.rosetta.core.select.jump_selector.InterchainJumpSelector())
        return mmf

    def get_mod_FastRelax(self, cycles: int = 1, weight: float = 10.0) -> pyrosetta.rosetta.protocols.moves.Mover:
        """
        This is not the usual fastRelax. It uses a modded minimiser protocol!
        No repacking.

        :param cycles: number of cycles
        :param weight: 10 is strict. 5 is decent. 1 is traditional.
        :return:
        """
        scorefxn = self._get_scorefxn("ref2015_cart")
        movemap = self._get_movemap()
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        v = pyrosetta.rosetta.std.vector_std_string(['repeat %%nrepeats%%',
                                                     f'coord_cst_weight {weight}',
                                                     'scale:fa_rep 0.092',
                                                     'min 0.01',
                                                     'scale:fa_rep 0.323',
                                                     'min 0.01',
                                                     'scale:fa_rep 0.633',
                                                     'min 0.01',
                                                     'scale:fa_rep 1',
                                                     'min 0.00001',
                                                     'accept_to_best',
                                                     'endrepeat'])
        relax.set_script_from_lines(v)
        relax.set_movemap(movemap)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.cartesian(True)
        relax.constrain_relax_to_start_coords(True)  # set native causes a segfault.
        return relax

    def get_old_FastRelax(self, cycles=1) -> pyrosetta.rosetta.protocols.moves.Mover:
        """

        :param cycles:
        :return:
        """
        scorefxn = self._get_scorefxn("ref2015_cart")
        movemap = self._get_movemap()
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.set_movemap(movemap)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.cartesian(True)
        relax.constrain_relax_to_start_coords(True)  # set native causes a segfault.
        return relax

    def get_MinMover(self) -> pyrosetta.rosetta.protocols.moves.Mover:
        """
        This gets stuck.

        :return: Move
        """
        scorefxn = self._get_scorefxn("ref2015_cart")
        movemap = self._get_movemap()
        minimizer = pyrosetta.rosetta.protocols.minimization_packing.MinMover(movemap,
                                                                              scorefxn,
                                                                              'lbfgs_armijo_nonmonotone',
                                                                              0.001,
                                                                              True)
        minimizer.cartesian(True)
        return minimizer

    def get_PertMinMover(self) -> pyrosetta.rosetta.protocols.moves.Mover:
        """
        This makes a mess of the structure.

        :return: Move
        """
        scorefxn = self._get_scorefxn(name="ref2015_cart")
        mmf = self._get_movemapfactory()
        minimizer = pyrosetta.rosetta.protocols.minimization_packing.PertMinMover()
        minimizer.movemap_factory(mmf)
        minimizer.scorefxn(scorefxn)
        minimizer.sc_only(False)
        minimizer.pert_size(1)
        minimizer.uniform(True)
        # minimizer.cartesian is True already thanks to movemap factory.
        return minimizer

    def repack_neighbors(self) -> None:
        # get score function. but relax the coordinate constraint.
        cc = self.coordinate_constraint
        self.coordinate_constraint = 0
        scorefxn = self._get_scorefxn("ref2015")
        self.coordinate_constraint = cc
        # get neighbourhood
        self._get_selector(ligand_only=True)
        NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
        ns = NeighborhoodResidueSelector(self._get_selector(ligand_only=True), distance=7,
                                         include_focus_in_subset=False)
        ## repack
        operation = pyrosetta.rosetta.core.pack.task.operation
        allow = operation.RestrictToRepackingRLT()
        restrict_to_focus = operation.OperateOnResidueSubset(allow, ns, True)
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        tf.push_back(operation.PreventRepacking())
        tf.push_back(operation.DisallowIfNonnative())
        tf.push_back(restrict_to_focus)
        packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
        packer.task_factory(tf)
        packer.apply(self.pose)

    def ligand_score(self):
        lig_pos = self.ligand_residue[0]
        # no constraints here
        scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("ref2015")
        scorefxn(self.pose)
        data = self.pose.energies().residue_total_energies_array()  # structured numpy array
        # this stupid line actually solves a race condition...
        assert data.shape[0] >= lig_pos - 1, f'Ligand {lig_pos} was lost from the pose? size={data.shape}'
        i = lig_pos - 1  ##pose numbering is fortran style. while python is C++
        sfxd = {data.dtype.names[j]: data[i][j] for j in range(len(data.dtype))}
        return {'MMFF_ligand': self.MMFF_score(delta=True),
                'holo_ref2015': scorefxn(self.pose),
                'ligand_ref2015': sfxd,
                **self.score_split()}

    def minimise(self, cycles: int = 10):
        #self.repack_neighbors()
        mover = self.get_mod_FastRelax(cycles)
        # mover = self.get_PertMinMover()
        # mover = self.get_MinMover()
        mover.apply(self.pose)

    def score_split(self, repack=False):
        split_pose = pyrosetta.Pose()
        split_pose.assign(self.pose)
        ResidueVector = pyrosetta.rosetta.core.select.residue_selector.ResidueVector
        x = self._get_selector(ligand_only=True).apply(split_pose)
        lig_pos = list(ResidueVector(self._get_selector(ligand_only=True).apply(split_pose)))[0]
        if self.pose.residue(lig_pos).connect_map_size() > 0:
            cys_pos = self.pose.residue(lig_pos).connect_map(1).resid()
            # RESCON: 305 LIG n-conn= 1 n-poly= 0 n-nonpoly= 1 conn# 1 22 145 3
            split_pose.conformation().sever_chemical_bond(seqpos1=cys_pos, res1_resconn_index=3, seqpos2=lig_pos,
                                                          res2_resconn_index=1)
        xyz = pyrosetta.rosetta.numeric.xyzVector_double_t()
        xyz.x = 500.0
        xyz.y = 0.0
        xyz.z = 0.0
        for a in range(1, split_pose.residue(lig_pos).natoms() + 1):
            split_pose.residue(lig_pos).set_xyz(a, split_pose.residue(lig_pos).xyz(a) + xyz)
        scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("ref2015")
        if repack:
            # get neighbourhood
            self._get_selector(ligand_only=True)
            NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
            ns = NeighborhoodResidueSelector(self._get_selector(ligand_only=True), distance=7,
                                             include_focus_in_subset=True)
            nsv = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(ns.apply(self.pose))
            fake_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(nsv)
            ## repack
            operation = pyrosetta.rosetta.core.pack.task.operation
            allow = operation.RestrictToRepackingRLT()
            restrict_to_focus = operation.OperateOnResidueSubset(allow, fake_selector, True)
            tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
            tf.push_back(operation.PreventRepacking())
            tf.push_back(restrict_to_focus)
            packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
            packer.task_factory(tf)
            packer.apply(split_pose)
        x = scorefxn(split_pose)
        b = scorefxn(self.pose)
        return {'xyz_unbound': x, 'xyz_bound': b, 'xyz_∆∆G': b - x}
