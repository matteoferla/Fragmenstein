"""File upload endpoints for template PDB and hit molecules."""

from fastapi import APIRouter, HTTPException, UploadFile

from ..database import get_db
from ..models.session import get_session, update_session
from ..schemas.upload import HitInfo, HitsResponse, UploadResponse
from ..services import file_manager
from ..services.molecule_service import mol_info, mol_to_molblock, parse_mol_file

router = APIRouter(prefix="/api/sessions/{session_id}", tags=["upload"])


@router.post("/template", response_model=UploadResponse)
async def upload_template(session_id: str, file: UploadFile):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    if not file.filename or not file.filename.lower().endswith(".pdb"):
        raise HTTPException(status_code=400, detail="Template must be a PDB file")

    content = await file.read()
    pdb_text = content.decode("utf-8", errors="replace")

    # Save file to disk and PDB text to database
    file_manager.save_upload(session_id, file.filename, content)
    update_session(
        session_id,
        template_filename=file.filename,
        template_pdb=pdb_text,
        status="template_uploaded",
    )

    return UploadResponse(filename=file.filename, message="Template uploaded successfully")


@router.post("/hits", response_model=UploadResponse)
async def upload_hits(session_id: str, files: list[UploadFile]):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    total_mols = 0
    for file in files:
        if not file.filename:
            continue

        content = await file.read()
        filepath = file_manager.save_upload(session_id, file.filename, content, subdir="hits")

        # Parse molecules from file
        try:
            mols = parse_mol_file(filepath)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))

        # Store each molecule in the database
        with get_db() as db:
            for mol in mols:
                info = mol_info(mol)
                mol_block = mol_to_molblock(mol)
                db.execute(
                    "INSERT INTO hits (session_id, name, filename, smiles, num_atoms, mol_block) "
                    "VALUES (?, ?, ?, ?, ?, ?)",
                    (session_id, info["name"], file.filename, info["smiles"], info["num_atoms"], mol_block),
                )
                total_mols += 1

    if total_mols > 0:
        update_session(session_id, status="hits_uploaded")

    return UploadResponse(
        filename=f"{len(files)} file(s)",
        message=f"Uploaded {total_mols} molecule(s) from {len(files)} file(s)",
    )


@router.get("/hits", response_model=HitsResponse)
def list_hits(session_id: str):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    with get_db() as db:
        rows = db.execute(
            "SELECT name, filename, smiles, num_atoms FROM hits WHERE session_id = ?",
            (session_id,),
        ).fetchall()

    hits = [
        HitInfo(
            name=r["name"],
            filename=r["filename"],
            smiles=r["smiles"],
            num_atoms=r["num_atoms"],
        )
        for r in rows
    ]
    return HitsResponse(hits=hits, count=len(hits))


@router.get("/template/pdb")
def get_template_pdb(session_id: str):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")
    if session.template_pdb is None:
        raise HTTPException(status_code=404, detail="No template uploaded")
    return {"pdb": session.template_pdb}
