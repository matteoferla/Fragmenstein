"""Serve MolBlock data for 3D viewer and 2D depictions."""

from fastapi import APIRouter, HTTPException
from fastapi.responses import Response
from rdkit import Chem
from rdkit.Chem import Draw

from ..database import get_db
from ..models.session import get_session

router = APIRouter(tags=["molecules"])


@router.get("/api/depict")
def depict_smiles(smiles: str, width: int = 250, height: int = 180):
    """Render a SMILES string as a PNG 2D structure depiction."""
    width = min(width, 1024)
    height = min(height, 1024)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES")
    svg = Draw.MolToImage(mol, size=(width, height))
    import io

    buf = io.BytesIO()
    svg.save(buf, format="PNG")
    buf.seek(0)
    return Response(
        content=buf.read(),
        media_type="image/png",
        headers={"Cache-Control": "public, max-age=86400"},
    )


@router.get("/api/sessions/{session_id}/hits/{hit_name}/mol")
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


@router.get("/api/sessions/{session_id}/hits/molblocks")
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

    return {"hits": [{"name": r["name"], "mol_block": r["mol_block"]} for r in rows]}
