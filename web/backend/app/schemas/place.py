"""Pydantic schemas for place endpoints."""

from pydantic import BaseModel


class PlaceRequest(BaseModel):
    victor_type: str = "Wictor"
    n_cores: int = -1
    timeout: int = 240
    merging_mode: str = "expansion"
    run_plip: bool = False
    covalent_resi: str | None = None
    source_job_id: str | None = None  # combine job to source hits from


class PlaceResultRow(BaseModel):
    index: int
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
    LE: float | None = None
    outcome: str | None = None
    const_ratio: float | None = None
    hit_names: list[str] | None = None
