"""File upload endpoints for template PDB and hit molecules."""

import asyncio
import logging

from fastapi import APIRouter, HTTPException, UploadFile

from ..database import get_db
from ..models.session import get_session, update_session
from ..schemas.upload import HitInfo, HitsResponse, TemplatePrepRequest, UploadResponse
from ..services import file_manager, job_manager
from ..services.molecule_service import mol_info, mol_to_molblock, parse_mol_file, parse_smi_file

log = logging.getLogger(__name__)

router = APIRouter(prefix="/api/sessions/{session_id}", tags=["upload"])


@router.post("/template", response_model=UploadResponse)
async def upload_template(session_id: str, file: UploadFile):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    if not file.filename or not file.filename.lower().endswith(".pdb"):
        raise HTTPException(status_code=400, detail="Template must be a PDB file")

    content = await file.read()
    pdb_text = content.decode("utf-8", errors="replace")

    # Save file to disk and PDB text to database
    file_manager.save_upload(session_id, file.filename, content)
    update_session(
        session_id,
        template_filename=file.filename,
        template_pdb=pdb_text,
        status="template_uploaded",
    )

    return UploadResponse(filename=file.filename, message="Template uploaded successfully")


@router.post("/hits", response_model=UploadResponse)
async def upload_hits(
    session_id: str,
    files: list[UploadFile],
    ligand_resn: str | None = None,
    proximity_bonding: bool = True,
):
    """Upload hit molecule files.

    For PDB files containing crystal structures, set ligand_resn to the
    3-letter residue name of the ligand (e.g. 'LIG', 'UNL') to extract
    just the ligand from the protein structure.
    """
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    # First pass: save all files and parse .smi for SMILES mapping
    smiles_map: dict[str, str] = {}
    saved_files: list[tuple[str, object]] = []

    for file in files:
        if not file.filename:
            continue
        content = await file.read()
        filepath = file_manager.save_upload(session_id, file.filename, content, subdir="hits")

        if file.filename.lower().endswith(".smi"):
            smiles_map.update(parse_smi_file(filepath))
        else:
            saved_files.append((file.filename, filepath))

    # Second pass: parse molecule files using the SMILES map for PDB bond order correction
    total_mols = 0
    for filename, filepath in saved_files:
        try:
            mols = parse_mol_file(
                filepath,
                ligand_resn=ligand_resn,
                proximity_bonding=proximity_bonding,
                smiles_map=smiles_map if smiles_map else None,
            )
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))

        # Store each molecule in the database
        with get_db() as db:
            for mol in mols:
                info = mol_info(mol)
                mol_block = mol_to_molblock(mol)
                db.execute(
                    "INSERT INTO hits (session_id, name, filename, smiles, num_atoms, mol_block) "
                    "VALUES (?, ?, ?, ?, ?, ?)",
                    (session_id, info["name"], filename, info["smiles"], info["num_atoms"], mol_block),
                )
                total_mols += 1

    if total_mols > 0:
        update_session(session_id, status="hits_uploaded")

    return UploadResponse(
        filename=f"{len(files)} file(s)",
        message=f"Uploaded {total_mols} molecule(s) from {len(files)} file(s)",
    )


@router.get("/hits", response_model=HitsResponse)
def list_hits(session_id: str):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    with get_db() as db:
        rows = db.execute(
            "SELECT name, filename, smiles, num_atoms FROM hits WHERE session_id = ?",
            (session_id,),
        ).fetchall()

    hits = [
        HitInfo(
            name=r["name"],
            filename=r["filename"],
            smiles=r["smiles"],
            num_atoms=r["num_atoms"],
        )
        for r in rows
    ]
    return HitsResponse(hits=hits, count=len(hits))


@router.post("/template/clean")
def clean_template(session_id: str, remove_residues: str = "HOH"):
    """Remove unwanted residues (waters, ions, etc.) from the template PDB.

    remove_residues: space/comma separated 3-letter codes, e.g. 'HOH SO4 CL'
    """
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    if session.template_pdb is None:
        raise HTTPException(status_code=400, detail="No template uploaded")

    import re
    resnames = [x for x in re.sub(r"[\W]", " ", remove_residues).split() if x]
    if not resnames:
        return {"message": "No residues to remove", "pdb_lines": len(session.template_pdb.split("\n"))}

    # Filter PDB lines — remove ATOM/HETATM lines matching unwanted residue names
    lines = session.template_pdb.split("\n")
    cleaned = []
    removed_count = 0
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if resname in resnames:
                removed_count += 1
                continue
        cleaned.append(line)

    cleaned_pdb = "\n".join(cleaned)
    update_session(session_id, template_pdb=cleaned_pdb)

    return {
        "message": f"Removed {removed_count} atoms from residues: {', '.join(resnames)}",
        "removed_count": removed_count,
    }


@router.get("/template/pdb")
def get_template_pdb(session_id: str):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    if session.template_pdb is None:
        raise HTTPException(status_code=404, detail="No template uploaded")
    return {"pdb": session.template_pdb}


async def _run_prepare_task(job_id: str, session_id: str, config: TemplatePrepRequest):
    """Background task that runs template preparation."""
    job_manager.mark_running(job_id)
    try:
        from ..services.template_prep_service import run_template_prepare

        summary = await asyncio.to_thread(
            run_template_prepare,
            session_id=session_id,
            parameterize=config.parameterize,
            minimize=config.minimize,
            center_resi=config.center_resi,
            center_chain=config.center_chain,
            neighborhood_radius=config.neighborhood_radius,
            cycles=config.cycles,
            remove_residues=config.remove_residues,
        )
        job_manager.mark_completed(job_id, result_path=summary)
    except Exception as e:
        log.exception(f"Template preparation failed for session {session_id}")
        job_manager.mark_failed(job_id, error=str(e))


@router.post("/template/prepare")
async def prepare_template(session_id: str, config: TemplatePrepRequest):
    """Prepare template using PyRosetta: parameterize, minimize, remove residues."""
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    if session.template_pdb is None:
        raise HTTPException(status_code=400, detail="Upload a template first")

    try:
        import pyrosetta  # noqa: F401
    except ImportError:
        raise HTTPException(status_code=400, detail="PyRosetta is not installed")

    job = job_manager.start_job(session_id, "template_prepare", config.model_dump())
    asyncio.create_task(_run_prepare_task(job.id, session_id, config))
    return {"job_id": job.id}
