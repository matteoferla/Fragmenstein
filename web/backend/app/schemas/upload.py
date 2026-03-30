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


class TemplatePrepRequest(BaseModel):
    parameterize: bool = True
    minimize: bool = False
    center_resi: int | None = None
    center_chain: str = "A"
    neighborhood_radius: float = 8.0
    cycles: int = 3
    remove_residues: list[str] = []
