"""Pydantic schemas for single-molecule Victor endpoints."""

from pydantic import BaseModel


class SingleCombineRequest(BaseModel):
    hit_names: list[str]
    victor_type: str = "Wictor"
    warhead_harmonisation: str = "first"
    joining_cutoff: float = 5.0
    covalent_resi: str | None = None


class SinglePlaceRequest(BaseModel):
    smiles: str
    hit_names: list[str]
    victor_type: str = "Wictor"
    merging_mode: str = "expansion"
    covalent_resi: str | None = None


class SingleResultResponse(BaseModel):
    name: str
    smiles: str | None = None
    error: str | None = None
    mode: str | None = None
    ddG: float | None = None
    dG_bound: float | None = None
    dG_unbound: float | None = None
    comRMSD: float | None = None
    N_constrained_atoms: int | None = None
    N_unconstrained_atoms: int | None = None
    runtime: float | None = None
    mol_block: str | None = None
    unmin_mol_block: str | None = None
