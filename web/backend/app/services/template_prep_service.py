"""Template preparation using PyRosetta — parameterization, energy minimization, residue removal."""

import logging
import re
from pathlib import Path

from . import file_manager
from .combine_service import load_session_template
from ..models.session import update_session

log = logging.getLogger(__name__)


def run_template_prepare(
    session_id: str,
    parameterize: bool = True,
    minimize: bool = False,
    center_resi: int | None = None,
    center_chain: str = "A",
    neighborhood_radius: float = 8.0,
    cycles: int = 3,
    remove_residues: list[str] | None = None,
) -> str:
    """Prepare the template PDB using PyRosetta.

    Steps:
    1. Load into PyRosetta (parameterize novel ligands if needed)
    2. Optionally energy-minimize around a target residue (FastRelax)
    3. Optionally remove unwanted residues (HOH, ions, etc.)

    Updates session.template_pdb with the prepared PDB and returns a summary.
    """
    import pyrosetta
    pyrosetta.distributed.maybe_init(
        extra_options="-no_optH false -ex1 -ex2 -mute all "
                      "-ignore_unrecognized_res true -load_PDB_components false "
                      "-ignore_waters true"
    )

    pdbblock = load_session_template(session_id)
    work_dir = file_manager.session_work_dir(session_id)
    steps_done = []

    # Step 1: Load pose (with parameterization of novel ligands)
    log.info(f"Loading template pose for session {session_id}")
    if parameterize:
        pose = _load_with_params(pdbblock, work_dir)
        steps_done.append("parameterized")
    else:
        pose = _pose_from_pdbstring(pdbblock)
        steps_done.append("loaded")

    # Step 2: Energy minimization
    if minimize and center_resi is not None:
        log.info(f"Running FastRelax: resi={center_resi}{center_chain}, radius={neighborhood_radius}, cycles={cycles}")
        _run_fast_relax(pose, center_resi, center_chain, neighborhood_radius, cycles)
        steps_done.append(f"minimized around {center_resi}{center_chain}")

    # Step 3: Remove residues
    if remove_residues:
        removed = _remove_residues(pose, remove_residues)
        if removed:
            steps_done.append(f"removed {', '.join(removed)}")

    # Extract prepared PDB
    prepared_pdb = _pose_to_pdbstring(pose)
    update_session(session_id, template_pdb=prepared_pdb)

    summary = "Template prepared: " + "; ".join(steps_done)
    log.info(summary)
    return summary


def _pose_from_pdbstring(pdbblock: str) -> "pyrosetta.Pose":
    """Create a Pose from a PDB string."""
    import pyrosetta
    from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring as _pdb_to_pose

    pose = pyrosetta.Pose()
    _pdb_to_pose(pose, pdbblock)
    return pose


def _load_with_params(pdbblock: str, work_dir: Path) -> "pyrosetta.Pose":
    """Load PDB into PyRosetta, generating params for novel ligands."""
    import pyrosetta

    # First try loading directly — works if all residues are known
    try:
        pose = _pose_from_pdbstring(pdbblock)
        return pose
    except RuntimeError:
        log.info("Direct loading failed, attempting parameterization of novel ligands")

    # Find HETATM residue names that aren't standard
    standard = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "ACE", "NME", "NH2"}
    hetatm_resnames = set()
    for line in pdbblock.splitlines():
        if line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname and resname not in standard:
                hetatm_resnames.add(resname)

    # Try generating params for each novel ligand using rdkit_to_params
    params_files = []
    try:
        from rdkit_to_params import Params
        for resname in hetatm_resnames:
            params_path = work_dir / f"{resname}.params"
            if params_path.exists():
                params_files.append(str(params_path))
                continue
            try:
                p = Params.from_smiles_w_pdbblock(pdbblock, resname)
                p.dump(str(params_path))
                params_files.append(str(params_path))
                log.info(f"Generated params for {resname}")
            except Exception as e:
                log.warning(f"Could not parameterize {resname}: {e}")
    except ImportError:
        log.warning("rdkit_to_params not installed, skipping parameterization")

    # Also check for any .params files already in the work directory
    for f in work_dir.glob("*.params"):
        path_str = str(f)
        if path_str not in params_files:
            params_files.append(path_str)

    from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring as _pdb_to_pose

    pose = pyrosetta.Pose()
    if params_files:
        rts = pose.conformation().modifiable_residue_type_set_for_conf()
        for pf in params_files:
            rts.read_files_for_base_residue_types(
                pyrosetta.rosetta.utility.vector1_string([pf])
            )
    new_pose = pyrosetta.Pose()
    _pdb_to_pose(new_pose, pdbblock)
    return new_pose


def _run_fast_relax(
    pose: "pyrosetta.Pose",
    center_resi: int,
    center_chain: str,
    neighborhood_radius: float,
    cycles: int,
):
    """Run FastRelax around a target residue."""
    import pyrosetta
    from pyrosetta.rosetta.core.select import residue_selector

    center_index = pose.pdb_info().pdb2pose(res=center_resi, chain=center_chain)
    if center_index == 0:
        raise ValueError(f"Residue {center_resi}{center_chain} not found in the pose")

    resi_sele = residue_selector.ResidueIndexSelector(center_index)
    neighbor_sele = residue_selector.NeighborhoodResidueSelector(
        resi_sele, distance=neighborhood_radius, include_focus_in_subset=True,
    )

    subset = neighbor_sele.apply(pose)
    movemap = pyrosetta.MoveMap()
    movemap.set_bb(subset)
    movemap.set_chi(subset)

    scorefxn = pyrosetta.get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    relax.set_movemap(movemap)
    relax.apply(pose)


def _remove_residues(pose: "pyrosetta.Pose", residue_names: list[str]) -> list[str]:
    """Remove residues by 3-letter name. Returns list of names actually removed."""
    from pyrosetta.rosetta.core.select import residue_selector, get_residues_from_subset

    removed = []
    for resn in residue_names:
        resn = resn.strip().upper()
        if not resn:
            continue
        sele = residue_selector.ResidueNameSelector()
        sele.set_residue_name3(resn)
        try:
            subset = sele.apply(pose)
        except RuntimeError:
            continue

        indices = list(get_residues_from_subset(subset))
        if not indices:
            continue

        for idx in sorted(indices, reverse=True):
            pose.delete_residue_slow(idx)
        removed.append(resn)

    return removed


def _pose_to_pdbstring(pose: "pyrosetta.Pose") -> str:
    """Convert a PyRosetta Pose back to a PDB string."""
    import pyrosetta
    buf = pyrosetta.rosetta.std.ostringstream()
    pose.dump_pdb(buf)
    return buf.str()
