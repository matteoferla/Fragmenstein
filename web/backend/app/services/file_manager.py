"""File storage and path resolution for sessions."""

import shutil
from pathlib import Path

from ..config import settings


def session_dir(session_id: str) -> Path:
    """Get the base directory for a session's files."""
    d = settings.upload_dir / session_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def session_work_dir(session_id: str) -> Path:
    """Get the work directory for a session (where Victor writes output)."""
    d = settings.work_dir / session_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def template_path(session_id: str) -> Path:
    """Path to the uploaded template PDB file."""
    return session_dir(session_id) / "template.pdb"


def hits_dir(session_id: str) -> Path:
    """Directory for uploaded hit files."""
    d = session_dir(session_id) / "hits"
    d.mkdir(parents=True, exist_ok=True)
    return d


def results_dir(session_id: str, job_type: str) -> Path:
    """Directory for job results."""
    d = session_work_dir(session_id) / job_type
    d.mkdir(parents=True, exist_ok=True)
    return d


def save_upload(session_id: str, filename: str, content: bytes, subdir: str = "") -> Path:
    """Save an uploaded file and return its path."""
    if subdir:
        dest = session_dir(session_id) / subdir
        dest.mkdir(parents=True, exist_ok=True)
    else:
        dest = session_dir(session_id)
    filepath = dest / filename
    filepath.write_bytes(content)
    return filepath


def cleanup_session(session_id: str):
    """Remove all files for a session."""
    for d in [session_dir(session_id), session_work_dir(session_id)]:
        if d.exists():
            shutil.rmtree(d)
