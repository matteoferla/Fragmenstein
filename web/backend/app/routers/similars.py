"""Similars endpoints — SmallWorld, manual SMILES, CSV/Excel upload, PubChem."""

import asyncio
import logging

from fastapi import APIRouter, HTTPException, UploadFile
from pydantic import BaseModel

from ..models.job import get_job
from ..models.session import get_session
from ..services import file_manager, job_manager
from ..models.job import get_jobs_for_session
from ..services.similars_service import (
    filter_library_against_mergers,
    parse_manual_smiles,
    parse_uploaded_file,
    run_pubchem_search,
    run_similars_search,
)

log = logging.getLogger(__name__)

router = APIRouter(prefix="/api/sessions/{session_id}", tags=["similars"])


# ── SmallWorld (existing) ─────────────────────────────────────────────

class SimilarsRequest(BaseModel):
    combine_job_id: str | None = None
    top_n: int = 100
    dist: int = 25
    length: int = 200
    db: str = "REAL_dataset"
    outcome_filter: str = "acceptable"


async def _run_similars_task(job_id: str, session_id: str, config: SimilarsRequest):
    job_manager.mark_running(job_id)
    try:
        # Resolve combine job — from request or find latest completed
        cjid = config.combine_job_id
        if not cjid:
            combine_jobs = [j for j in get_jobs_for_session(session_id)
                            if j.type == "combine" and j.status == "completed" and j.result_path]
            if not combine_jobs:
                raise ValueError("No completed combine job found for this session")
            cjid = combine_jobs[0].id
        combine_job = get_job(cjid)
        if combine_job is None or combine_job.result_path is None:
            raise ValueError("Combine job not found or not completed")
        result_path = await asyncio.to_thread(
            run_similars_search,
            session_id=session_id,
            combine_result_path=combine_job.result_path,
            top_n=config.top_n,
            dist=config.dist,
            length=config.length,
            db=config.db,
            outcome_filter=config.outcome_filter,
        )
        job_manager.mark_completed(job_id, result_path=str(result_path))
    except Exception as e:
        log.exception(f"SmallWorld search failed for session {session_id}")
        job_manager.mark_failed(job_id, error=f"{type(e).__name__}: {e}")


@router.post("/similars")
async def start_similars(session_id: str, config: SimilarsRequest):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    job = job_manager.start_job(session_id, "similars", config.model_dump())
    asyncio.create_task(_run_similars_task(job.id, session_id, config))
    return {"job_id": job.id}


# ── Manual SMILES paste ──────────────────────────────────────────────

class ManualSmilesRequest(BaseModel):
    smiles_text: str


@router.post("/similars/manual")
def manual_smiles(session_id: str, req: ManualSmilesRequest):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    if not req.smiles_text.strip():
        raise HTTPException(status_code=400, detail="No SMILES provided")

    result_path, valid, invalid = parse_manual_smiles(session_id, req.smiles_text)

    # Create a completed job record so Place can reference it
    job = job_manager.start_job(session_id, "similars", {"source": "manual"})
    job_manager.mark_completed(job.id, result_path=str(result_path))

    return {
        "job_id": job.id,
        "valid": valid,
        "invalid": invalid,
        "message": f"Loaded {valid} compound(s)" + (f" ({invalid} invalid SMILES skipped)" if invalid else ""),
    }


# ── CSV/Excel upload ─────────────────────────────────────────────────

@router.post("/similars/upload")
async def upload_similars(session_id: str, file: UploadFile):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    if not file.filename:
        raise HTTPException(status_code=400, detail="No file provided")

    content = await file.read()
    filepath = file_manager.save_upload(session_id, file.filename, content, subdir="similars")

    try:
        result_path, valid, invalid = parse_uploaded_file(session_id, filepath)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    job = job_manager.start_job(session_id, "similars", {"source": "upload", "filename": file.filename})
    job_manager.mark_completed(job.id, result_path=str(result_path))

    return {
        "job_id": job.id,
        "valid": valid,
        "invalid": invalid,
        "message": f"Loaded {valid} compound(s) from {file.filename}" + (f" ({invalid} invalid SMILES skipped)" if invalid else ""),
    }


# ── Upload + filter against mergers ──────────────────────────────────

@router.post("/similars/upload-filter")
async def upload_and_filter(
    session_id: str,
    file: UploadFile,
    top_n: int = 200,
    outcome_filter: str = "acceptable",
):
    """Upload a large compound library and filter by Tanimoto similarity to mergers.

    Uses Morgan fingerprints (radius=2) to find the top_n most similar compounds
    from the uploaded library relative to the combine step's mergers.
    """
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    if not file.filename:
        raise HTTPException(status_code=400, detail="No file provided")

    # Find completed combine job
    combine_jobs = [j for j in get_jobs_for_session(session_id)
                    if j.type == "combine" and j.status == "completed" and j.result_path]
    if not combine_jobs:
        raise HTTPException(status_code=400, detail="No completed combine job found. Run Combine first.")

    content = await file.read()
    filepath = file_manager.save_upload(session_id, file.filename, content, subdir="similars")

    # First parse to validate SMILES
    try:
        library_path, valid, invalid = parse_uploaded_file(session_id, filepath)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    if valid == 0:
        raise HTTPException(status_code=400, detail="No valid SMILES found in file")

    # Filter against mergers
    try:
        result_path, library_size, filtered_count, merger_count = await asyncio.to_thread(
            filter_library_against_mergers,
            session_id=session_id,
            library_path=library_path,
            combine_result_path=combine_jobs[0].result_path,
            top_n=top_n,
            outcome_filter=outcome_filter,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    job = job_manager.start_job(session_id, "similars", {
        "source": "upload-filter",
        "filename": file.filename,
        "library_size": library_size,
        "top_n": top_n,
    })
    job_manager.mark_completed(job.id, result_path=str(result_path))

    return {
        "job_id": job.id,
        "library_size": library_size,
        "filtered": filtered_count,
        "mergers_used": merger_count,
        "invalid": invalid,
        "message": f"Filtered {library_size} compounds → top {filtered_count} by Tanimoto similarity to {merger_count} mergers" +
                   (f" ({invalid} invalid SMILES skipped)" if invalid else ""),
    }


# ── PubChem similarity search ────────────────────────────────────────

class PubChemRequest(BaseModel):
    combine_job_id: str | None = None
    top_n: int = 50
    threshold: int = 80
    max_per_query: int = 20
    outcome_filter: str = "acceptable"


async def _run_pubchem_task(job_id: str, session_id: str, config: PubChemRequest):
    job_manager.mark_running(job_id)
    try:
        # Resolve combine job — from request or find latest completed
        cjid = config.combine_job_id
        if not cjid:
            combine_jobs = [j for j in get_jobs_for_session(session_id)
                            if j.type == "combine" and j.status == "completed" and j.result_path]
            if not combine_jobs:
                raise ValueError("No completed combine job found for this session")
            cjid = combine_jobs[0].id
        combine_job = get_job(cjid)
        if combine_job is None or combine_job.result_path is None:
            raise ValueError("Combine job not found or not completed")
        result_path = await asyncio.to_thread(
            run_pubchem_search,
            session_id=session_id,
            combine_result_path=combine_job.result_path,
            top_n=config.top_n,
            threshold=config.threshold,
            max_per_query=config.max_per_query,
            outcome_filter=config.outcome_filter,
        )
        job_manager.mark_completed(job_id, result_path=str(result_path))
    except Exception as e:
        log.exception(f"PubChem search failed for session {session_id}")
        job_manager.mark_failed(job_id, error=f"{type(e).__name__}: {e}")


@router.post("/similars/pubchem")
async def start_pubchem(session_id: str, config: PubChemRequest):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    job = job_manager.start_job(session_id, "similars", config.model_dump())
    asyncio.create_task(_run_pubchem_task(job.id, session_id, config))
    return {"job_id": job.id}
