"""Pydantic schemas for upload endpoints."""

from pydantic import BaseModel


class FileInfo(BaseModel):
    filename: str
    size: int


class UploadResponse(BaseModel):
    filename: str
    message: str


class HitInfo(BaseModel):
    name: str
    filename: str
    smiles: str | None = None
    num_atoms: int | None = None


class HitsResponse(BaseModel):
    hits: list[HitInfo]
    count: int
