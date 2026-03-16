"""Pydantic schemas for combine endpoints."""

from pydantic import BaseModel, Field


class CombineRequest(BaseModel):
    victor_type: str = "Wictor"
    n_cores: int = -1
    timeout: int = 240
    combination_size: int = 2
    permute: bool = True
    joining_cutoff: float = 5.0
    quick_reanimation: bool = False
    warhead_harmonisation: str = "first"
    run_plip: bool = False
    covalent_resi: str | None = None
    hit_names: list[str] | None = None  # subset of hits to combine; None = all


class CombineResultRow(BaseModel):
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
    percent_hybrid: float | None = None
    largest_ring: int | None = None
    N_HA: int | None = None
    N_rotatable_bonds: int | None = None
    hit_names: list[str] | None = None
