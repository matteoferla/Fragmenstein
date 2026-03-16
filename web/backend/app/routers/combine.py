"""Combine endpoint — starts Laboratory.combine() as a background job."""

import asyncio
import logging

from fastapi import APIRouter, HTTPException

from ..models.session import get_session
from ..schemas.combine import CombineRequest, CombineResultRow
from ..services import job_manager
from ..services.combine_service import run_combine
from ..services.result_serializer import dataframe_to_rows, load_dataframe

log = logging.getLogger(__name__)

router = APIRouter(prefix="/api/sessions/{session_id}", tags=["combine"])


async def _run_combine_task(job_id: str, session_id: str, config: CombineRequest):
    """Background task that runs combine and updates job status."""
    job_manager.mark_running(job_id)
    try:
        result_path = await asyncio.to_thread(
            run_combine,
            session_id=session_id,
            victor_type=config.victor_type,
            n_cores=config.n_cores,
            timeout=config.timeout,
            combination_size=config.combination_size,
            permute=config.permute,
            joining_cutoff=config.joining_cutoff,
            quick_reanimation=config.quick_reanimation,
            warhead_harmonisation=config.warhead_harmonisation,
            run_plip=config.run_plip,
            covalent_resi=config.covalent_resi,
            hit_names=config.hit_names,
        )
        job_manager.mark_completed(job_id, result_path=str(result_path))
    except Exception as e:
        log.exception(f"Combine failed for session {session_id}")
        job_manager.mark_failed(job_id, error=str(e))


@router.post("/combine")
async def start_combine(session_id: str, config: CombineRequest):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    if session.template_pdb is None:
        raise HTTPException(status_code=400, detail="Upload a template first")

    job = job_manager.start_job(session_id, "combine", config.model_dump())
    asyncio.create_task(_run_combine_task(job.id, session_id, config))
    return {"job_id": job.id}
