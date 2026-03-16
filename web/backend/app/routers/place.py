"""Place endpoint — starts Laboratory.place() as a background job."""

import asyncio
import logging
from pathlib import Path

from fastapi import APIRouter, HTTPException

from ..models.job import get_job
from ..models.session import get_session
from ..schemas.place import PlaceRequest
from ..services import job_manager
from ..services.place_service import build_place_queries_from_similars, run_place
from ..services.result_serializer import load_dataframe

log = logging.getLogger(__name__)

router = APIRouter(prefix="/api/sessions/{session_id}", tags=["place"])


async def _run_place_task(job_id: str, session_id: str, config: PlaceRequest):
    job_manager.mark_running(job_id)
    try:
        # Build queries from similars results
        if config.source_job_id:
            source_job = get_job(config.source_job_id)
            if source_job is None or source_job.result_path is None:
                raise ValueError("Source job not found or not completed")
            similars_df = load_dataframe(Path(source_job.result_path))
            queries = build_place_queries_from_similars(session_id, similars_df)
        else:
            raise ValueError("source_job_id is required for placement")

        result_path = await asyncio.to_thread(
            run_place,
            session_id=session_id,
            queries=queries,
            victor_type=config.victor_type,
            n_cores=config.n_cores,
            timeout=config.timeout,
            covalent_resi=config.covalent_resi,
        )
        job_manager.mark_completed(job_id, result_path=str(result_path))
    except Exception as e:
        log.exception(f"Place failed for session {session_id}")
        job_manager.mark_failed(job_id, error=str(e))


@router.post("/place")
async def start_place(session_id: str, config: PlaceRequest):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    job = job_manager.start_job(session_id, "place", config.model_dump())
    asyncio.create_task(_run_place_task(job.id, session_id, config))
    return {"job_id": job.id}
