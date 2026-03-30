"""Single-molecule Victor endpoints — run one combine or place."""

import asyncio
import logging

from fastapi import APIRouter, HTTPException

from ..models.session import get_session
from ..schemas.single_victor import SingleCombineRequest, SinglePlaceRequest, SingleResultResponse
from ..services.single_victor_service import run_single_combine, run_single_place

log = logging.getLogger(__name__)

router = APIRouter(prefix="/api/sessions/{session_id}/single", tags=["single"])


@router.post("/combine", response_model=SingleResultResponse)
async def single_combine(session_id: str, config: SingleCombineRequest):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    if session.template_pdb is None:
        raise HTTPException(status_code=400, detail="Upload a template first")

    try:
        result = await asyncio.to_thread(
            run_single_combine,
            session_id=session_id,
            hit_names=config.hit_names,
            victor_type=config.victor_type,
            warhead_harmonisation=config.warhead_harmonisation,
            joining_cutoff=config.joining_cutoff,
            covalent_resi=config.covalent_resi,
        )
        return SingleResultResponse(**result)
    except Exception as e:
        log.exception(f"Single combine failed for session {session_id}")
        return SingleResultResponse(name="error", error=str(e))


@router.post("/place", response_model=SingleResultResponse)
async def single_place(session_id: str, config: SinglePlaceRequest):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    if session.template_pdb is None:
        raise HTTPException(status_code=400, detail="Upload a template first")

    try:
        result = await asyncio.to_thread(
            run_single_place,
            session_id=session_id,
            smiles=config.smiles,
            hit_names=config.hit_names,
            victor_type=config.victor_type,
            merging_mode=config.merging_mode,
            covalent_resi=config.covalent_resi,
        )
        return SingleResultResponse(**result)
    except Exception as e:
        log.exception(f"Single place failed for session {session_id}")
        return SingleResultResponse(name="error", error=str(e))
