"""Job status and SSE progress stream endpoints."""

import io
from pathlib import Path

from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse

from ..models.job import get_job, get_jobs_for_session
from ..schemas.job import JobStatusResponse
from ..services import job_manager
from ..services.result_serializer import dataframe_to_rows, get_mol_block, load_dataframe, similars_dataframe_to_rows

router = APIRouter(prefix="/api/jobs", tags=["jobs"])


@router.get("/{job_id}", response_model=JobStatusResponse)
def get_status(job_id: str):
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return JobStatusResponse(
        id=job.id,
        session_id=job.session_id,
        type=job.type,
        status=job.status,
        progress=job.progress,
        message=job.message,
        error=job.error,
        created_at=job.created_at,
        started_at=job.started_at,
        completed_at=job.completed_at,
    )


@router.post("/{job_id}/cancel")
def cancel_job(job_id: str):
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.status in ("completed", "failed", "cancelled"):
        return {"status": job.status, "message": "Job already finished"}
    job_manager.mark_cancelled(job_id)
    return {"status": "cancelled", "message": "Job cancelled"}


@router.get("/{job_id}/stream")
async def stream_progress(job_id: str):
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return StreamingResponse(
        job_manager.subscribe(job_id),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no",
        },
    )


@router.get("/{job_id}/results")
def get_results(job_id: str):
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.status != "completed" or job.result_path is None:
        raise HTTPException(status_code=400, detail="Job not completed")

    df = load_dataframe(Path(job.result_path))
    if job.type == "similars":
        rows = similars_dataframe_to_rows(df)
    else:
        rows = dataframe_to_rows(df)
    return {"results": rows, "count": len(rows)}


@router.get("/{job_id}/results/{idx}/mol")
def get_result_mol(job_id: str, idx: int, mol_type: str = "minimized"):
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.result_path is None:
        raise HTTPException(status_code=400, detail="No results available")

    df = load_dataframe(Path(job.result_path))
    mol_block = get_mol_block(df, idx, mol_type)
    if mol_block is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {"mol_block": mol_block}


@router.get("/{job_id}/results/download")
def download_results(job_id: str, format: str = "csv"):
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.result_path is None:
        raise HTTPException(status_code=400, detail="No results available")

    df = load_dataframe(Path(job.result_path))

    if format == "csv":
        # Drop Mol columns for CSV export
        scalar_cols = [c for c in df.columns if not any(
            kw in str(c).lower() for kw in ["mol", "hit_mols", "binary"]
        ) and not isinstance(c, tuple)]
        buf = io.StringIO()
        df[scalar_cols].to_csv(buf, index=False)
        return StreamingResponse(
            iter([buf.getvalue()]),
            media_type="text/csv",
            headers={"Content-Disposition": f"attachment; filename={job.type}_results.csv"},
        )
    elif format == "sdf":
        from rdkit import Chem
        buf = io.StringIO()
        mol_col = None
        for col in ["minimized_mol", "unminimized_mol"]:
            if col in df.columns:
                mol_col = col
                break
        if mol_col is None:
            raise HTTPException(status_code=400, detail="No molecule column for SDF export")
        writer = Chem.SDWriter(buf)
        for idx, row in df.iterrows():
            mol = row.get(mol_col)
            if mol is not None and isinstance(mol, Chem.Mol):
                mol.SetProp("_Name", str(row.get("name", f"mol_{idx}")))
                if "∆∆G" in row.index:
                    val = row["∆∆G"]
                    if val is not None and str(val) != "nan":
                        mol.SetProp("ddG", f"{val:.3f}")
                if "outcome" in row.index and row["outcome"]:
                    mol.SetProp("outcome", str(row["outcome"]))
                writer.write(mol)
        writer.close()
        return StreamingResponse(
            iter([buf.getvalue()]),
            media_type="chemical/x-mdl-sdfile",
            headers={"Content-Disposition": f"attachment; filename={job.type}_results.sdf"},
        )
    else:
        raise HTTPException(status_code=400, detail=f"Unsupported format: {format}")
