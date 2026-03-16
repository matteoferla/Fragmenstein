"""SmallWorld search endpoint."""

import asyncio
import logging

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from ..models.job import get_job
from ..models.session import get_session
from ..services import job_manager
from ..services.similars_service import run_similars_search

log = logging.getLogger(__name__)

router = APIRouter(prefix="/api/sessions/{session_id}", tags=["similars"])


class SimilarsRequest(BaseModel):
    combine_job_id: str
    top_n: int = 100
    dist: int = 25
    length: int = 200
    db: str = "REAL_dataset"
    outcome_filter: str = "acceptable"


async def _run_similars_task(job_id: str, session_id: str, config: SimilarsRequest):
    job_manager.mark_running(job_id)
    try:
        combine_job = get_job(config.combine_job_id)
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
        job_manager.mark_failed(job_id, error=str(e))


@router.post("/similars")
async def start_similars(session_id: str, config: SimilarsRequest):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    job = job_manager.start_job(session_id, "similars", config.model_dump())
    asyncio.create_task(_run_similars_task(job.id, session_id, config))
    return {"job_id": job.id}
