"""Demo dataset endpoint — load MPro or Mac1 data with one click."""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from rdkit import Chem

from ..database import get_db
from ..models.session import get_session, update_session
from ..services.molecule_service import mol_info, mol_to_molblock

router = APIRouter(prefix="/api/sessions/{session_id}", tags=["demo"])


class DemoResponse(BaseModel):
    template_filename: str
    hit_count: int
    message: str


@router.post("/demo/{dataset}", response_model=DemoResponse)
def load_demo(session_id: str, dataset: str, n_hits: int = 5):
    session = get_session(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    from fragmenstein.demo import MPro, Mac1

    datasets = {"mpro": MPro, "mac1": Mac1}
    ds = datasets.get(dataset.lower())
    if ds is None:
        raise HTTPException(status_code=400, detail=f"Unknown dataset: {dataset}. Use 'mpro' or 'mac1'.")

    # Load template
    template_pdb = ds.get_template()
    template_name = f"{dataset}_template.pdb"
    update_session(session_id, template_filename=template_name, template_pdb=template_pdb, status="template_uploaded")

    # Load hits
    hits = ds.get_n_filtered_mols(n_hits)
    with get_db() as db:
        # Clear existing hits for this session
        db.execute("DELETE FROM hits WHERE session_id = ?", (session_id,))
        for mol in hits:
            info = mol_info(mol)
            mol_block = mol_to_molblock(mol)
            db.execute(
                "INSERT INTO hits (session_id, name, filename, smiles, num_atoms, mol_block) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (session_id, info["name"], f"{dataset}_demo.sdf", info["smiles"], info["num_atoms"], mol_block),
            )

    update_session(session_id, status="hits_uploaded")

    return DemoResponse(
        template_filename=template_name,
        hit_count=len(hits),
        message=f"Loaded {dataset.upper()} demo: template + {len(hits)} hits",
    )
