"""Monster standalone endpoints — chemistry-only combine/place without protein."""

import asyncio
import logging

from fastapi import APIRouter, HTTPException

from ..models.session import get_session
from ..schemas.monster import MonsterCombineRequest, MonsterPlaceRequest, MonsterResultResponse
from ..services.monster_service import run_monster_combine, run_monster_place

log = logging.getLogger(__name__)

router = APIRouter(prefix="/api/sessions/{session_id}/monster", tags=["monster"])


@router.post("/combine", response_model=MonsterResultResponse)
async def monster_combine(session_id: str, config: MonsterCombineRequest):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    try:
        result = await asyncio.to_thread(
            run_monster_combine,
            session_id=session_id,
            joining_cutoff=config.joining_cutoff,
            hit_names=config.hit_names,
        )
        return MonsterResultResponse(**result)
    except Exception as e:
        log.exception(f"Monster combine failed for session {session_id}")
        return MonsterResultResponse(name="error", error=str(e))


@router.post("/place", response_model=MonsterResultResponse)
async def monster_place(session_id: str, config: MonsterPlaceRequest):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    try:
        result = await asyncio.to_thread(
            run_monster_place,
            session_id=session_id,
            smiles=config.smiles,
            hit_names=config.hit_names,
        )
        return MonsterResultResponse(**result)
    except Exception as e:
        log.exception(f"Monster place failed for session {session_id}")
        return MonsterResultResponse(name="error", error=str(e))
