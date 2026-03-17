"""Pydantic schemas for Monster standalone endpoints."""

from pydantic import BaseModel


class MonsterCombineRequest(BaseModel):
    joining_cutoff: float = 5.0
    hit_names: list[str] | None = None


class MonsterPlaceRequest(BaseModel):
    smiles: str
    hit_names: list[str] | None = None


class MonsterResultResponse(BaseModel):
    name: str
    smiles: str | None = None
    num_atoms: int | None = None
    mol_block: str | None = None
    error: str | None = None
