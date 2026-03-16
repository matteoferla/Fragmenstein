"""Serve MolBlock data for 3D viewer."""

from fastapi import APIRouter, HTTPException

from ..database import get_db
from ..models.session import get_session

router = APIRouter(prefix="/api/sessions/{session_id}", tags=["molecules"])


@router.get("/hits/{hit_name}/mol")
def get_hit_molblock(session_id: str, hit_name: str):
    """Get MolBlock for a specific hit molecule."""
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    with get_db() as db:
        row = db.execute(
            "SELECT mol_block FROM hits WHERE session_id = ? AND name = ?",
            (session_id, hit_name),
        ).fetchone()

    if row is None:
        raise HTTPException(status_code=404, detail="Hit not found")
    return {"mol_block": row["mol_block"]}


@router.get("/hits/molblocks")
def get_all_hit_molblocks(session_id: str):
    """Get MolBlocks for all hit molecules in a session."""
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    with get_db() as db:
        rows = db.execute(
            "SELECT name, mol_block FROM hits WHERE session_id = ?",
            (session_id,),
        ).fetchall()

    return {
        "hits": [{"name": r["name"], "mol_block": r["mol_block"]} for r in rows]
    }
