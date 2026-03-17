"""Pydantic schemas for session endpoints."""

from pydantic import BaseModel
from typing import Optional


class CreateSessionRequest(BaseModel):
    name: str = ""


class SessionResponse(BaseModel):
    id: str
    created_at: str
    name: str
    status: str
    template_filename: Optional[str] = None
    hit_count: int = 0
