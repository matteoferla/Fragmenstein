"""Session model — CRUD operations on the sessions table."""

import json
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional

from ..database import get_db


@dataclass
class Session:
    id: str
    created_at: str
    name: str
    status: str
    config_json: str = "{}"
    template_filename: Optional[str] = None
    template_pdb: Optional[str] = None

    @property
    def config(self) -> dict:
        return json.loads(self.config_json)


def create_session(name: str = "") -> Session:
    sid = str(uuid.uuid4())
    with get_db() as db:
        db.execute(
            "INSERT INTO sessions (id, name) VALUES (?, ?)",
            (sid, name),
        )
        row = db.execute("SELECT * FROM sessions WHERE id = ?", (sid,)).fetchone()
    return Session(**dict(row))


def get_session(session_id: str) -> Optional[Session]:
    with get_db() as db:
        row = db.execute("SELECT * FROM sessions WHERE id = ?", (session_id,)).fetchone()
    if row is None:
        return None
    return Session(**dict(row))


def list_sessions() -> list[Session]:
    with get_db() as db:
        rows = db.execute(
            "SELECT * FROM sessions ORDER BY created_at DESC"
        ).fetchall()
    return [Session(**dict(r)) for r in rows]


def update_session(session_id: str, **kwargs) -> Optional[Session]:
    allowed = {"name", "status", "config_json", "template_filename", "template_pdb"}
    updates = {k: v for k, v in kwargs.items() if k in allowed}
    if not updates:
        return get_session(session_id)
    set_clause = ", ".join(f"{k} = ?" for k in updates)
    values = list(updates.values()) + [session_id]
    with get_db() as db:
        db.execute(f"UPDATE sessions SET {set_clause} WHERE id = ?", values)
    return get_session(session_id)


def delete_session(session_id: str) -> bool:
    with get_db() as db:
        cursor = db.execute("DELETE FROM sessions WHERE id = ?", (session_id,))
    return cursor.rowcount > 0
