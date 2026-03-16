"""Job model — CRUD operations on the jobs table."""

import json
import uuid
from dataclasses import dataclass
from datetime import datetime
from typing import Optional

from ..database import get_db


@dataclass
class Job:
    id: str
    session_id: str
    type: str
    status: str
    progress: float
    message: str
    error: Optional[str]
    config_json: str
    result_path: Optional[str]
    created_at: str
    started_at: Optional[str]
    completed_at: Optional[str]

    @property
    def config(self) -> dict:
        return json.loads(self.config_json)


def create_job(session_id: str, job_type: str, config: dict | None = None) -> Job:
    jid = str(uuid.uuid4())
    config_json = json.dumps(config or {})
    with get_db() as db:
        db.execute(
            "INSERT INTO jobs (id, session_id, type, config_json) VALUES (?, ?, ?, ?)",
            (jid, session_id, job_type, config_json),
        )
        row = db.execute("SELECT * FROM jobs WHERE id = ?", (jid,)).fetchone()
    return Job(**dict(row))


def get_job(job_id: str) -> Optional[Job]:
    with get_db() as db:
        row = db.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
    if row is None:
        return None
    return Job(**dict(row))


def get_jobs_for_session(session_id: str) -> list[Job]:
    with get_db() as db:
        rows = db.execute(
            "SELECT * FROM jobs WHERE session_id = ? ORDER BY created_at DESC",
            (session_id,),
        ).fetchall()
    return [Job(**dict(r)) for r in rows]


def update_job(job_id: str, **kwargs) -> Optional[Job]:
    allowed = {"status", "progress", "message", "error", "result_path", "started_at", "completed_at"}
    updates = {k: v for k, v in kwargs.items() if k in allowed}
    if not updates:
        return get_job(job_id)
    set_clause = ", ".join(f"{k} = ?" for k in updates)
    values = list(updates.values()) + [job_id]
    with get_db() as db:
        db.execute(f"UPDATE jobs SET {set_clause} WHERE id = ?", values)
    return get_job(job_id)
