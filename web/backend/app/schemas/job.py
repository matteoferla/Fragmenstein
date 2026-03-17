"""Pydantic schemas for job endpoints."""

from pydantic import BaseModel
from typing import Optional


class JobStatusResponse(BaseModel):
    id: str
    session_id: str
    type: str
    status: str
    progress: float
    message: str
    error: Optional[str] = None
    created_at: str
    started_at: Optional[str] = None
    completed_at: Optional[str] = None


class JobProgressEvent(BaseModel):
    progress: float
    status: str
    message: str
