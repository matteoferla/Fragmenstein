"""Session CRUD endpoints."""

from fastapi import APIRouter, HTTPException

from ..database import get_db
from ..models.job import get_jobs_for_session
from ..models.session import (
    create_session,
    delete_session,
    get_session,
    list_sessions,
)
from ..schemas.job import JobStatusResponse
from ..schemas.session import CreateSessionRequest, SessionResponse

router = APIRouter(prefix="/api/sessions", tags=["sessions"])


def _session_response(session) -> SessionResponse:
    with get_db() as db:
        hit_count = db.execute(
            "SELECT COUNT(*) as cnt FROM hits WHERE session_id = ?",
            (session.id,),
        ).fetchone()["cnt"]
    return SessionResponse(
        id=session.id,
        created_at=session.created_at,
        name=session.name,
        status=session.status,
        template_filename=session.template_filename,
        hit_count=hit_count,
    )


@router.post("", response_model=SessionResponse, status_code=201)
def create(req: CreateSessionRequest):
    session = create_session(name=req.name)
    return _session_response(session)


@router.get("", response_model=list[SessionResponse])
def list_all():
    sessions = list_sessions()
    return [_session_response(s) for s in sessions]


@router.get("/{session_id}", response_model=SessionResponse)
def get(session_id: str):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    return _session_response(session)


@router.get("/{session_id}/jobs", response_model=list[JobStatusResponse])
def list_jobs(session_id: str):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    jobs = get_jobs_for_session(session_id)
    return [
        JobStatusResponse(
            id=j.id, session_id=j.session_id, type=j.type, status=j.status,
            progress=j.progress, message=j.message, error=j.error,
            created_at=j.created_at, started_at=j.started_at, completed_at=j.completed_at,
        )
        for j in jobs
    ]


@router.delete("/{session_id}", status_code=204)
def delete(session_id: str):
    from ..services.file_manager import cleanup_session

    if not delete_session(session_id):
        raise HTTPException(status_code=404, detail="Session not found")
    cleanup_session(session_id)
