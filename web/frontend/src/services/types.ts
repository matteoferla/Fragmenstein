/** TypeScript interfaces for the Fragmenstein API. */

export interface SessionResponse {
  id: string;
  created_at: string;
  name: string;
  status: string;
  template_filename: string | null;
  hit_count: number;
}

export interface HitInfo {
  name: string;
  filename: string;
  smiles: string | null;
  num_atoms: number | null;
}

export interface HitsResponse {
  hits: HitInfo[];
  count: number;
}

export interface CombineRequest {
  victor_type?: string;
  n_cores?: number;
  timeout?: number;
  combination_size?: number;
  permute?: boolean;
  joining_cutoff?: number;
  quick_reanimation?: boolean;
  warhead_harmonisation?: string;
  run_plip?: boolean;
  covalent_resi?: string | null;
  hit_names?: string[] | null;
}

export interface PlaceRequest {
  victor_type?: string;
  n_cores?: number;
  timeout?: number;
  merging_mode?: string;
  run_plip?: boolean;
  use_originals?: boolean;
  covalent_resi?: string | null;
  source_job_id?: string | null;
}

export interface SimilarsRequest {
  combine_job_id: string;
  top_n?: number;
  dist?: number;
  length?: number;
  db?: string;
  outcome_filter?: string;
}

export interface ResultRow {
  index: number;
  name: string;
  smiles: string | null;
  error: string | null;
  mode: string | null;
  ddG: number | null;
  dG_bound: number | null;
  dG_unbound: number | null;
  comRMSD: number | null;
  N_constrained_atoms: number | null;
  N_unconstrained_atoms: number | null;
  runtime: number | null;
  LE: number | null;
  outcome: string | null;
  percent_hybrid: number | null;
  largest_ring: number | null;
  N_HA: number | null;
  N_rotatable_bonds: number | null;
  const_ratio: number | null;
  hit_names: string[] | null;
}

export interface ResultsResponse {
  results: ResultRow[];
  count: number;
}

export interface JobStatus {
  id: string;
  session_id: string;
  type: string;
  status: "pending" | "running" | "completed" | "failed";
  progress: number;
  message: string;
  error: string | null;
  created_at: string;
  started_at: string | null;
  completed_at: string | null;
}

export interface JobProgressEvent {
  progress: number;
  status: string;
  message: string;
}

export interface MonsterCombineRequest {
  joining_cutoff?: number;
  hit_names?: string[] | null;
}

export interface MonsterPlaceRequest {
  smiles: string;
  hit_names?: string[] | null;
}

export interface MonsterResult {
  name: string;
  smiles: string | null;
  num_atoms: number | null;
  mol_block: string | null;
  error: string | null;
}

export interface SingleCombineRequest {
  hit_names: string[];
  victor_type?: string;
  warhead_harmonisation?: string;
  joining_cutoff?: number;
  covalent_resi?: string | null;
}

export interface SinglePlaceRequest {
  smiles: string;
  hit_names: string[];
  victor_type?: string;
  merging_mode?: string;
  covalent_resi?: string | null;
}

export interface SingleResult {
  name: string;
  smiles: string | null;
  error: string | null;
  mode: string | null;
  ddG: number | null;
  dG_bound: number | null;
  dG_unbound: number | null;
  comRMSD: number | null;
  N_constrained_atoms: number | null;
  N_unconstrained_atoms: number | null;
  runtime: number | null;
  mol_block: string | null;
  unmin_mol_block: string | null;
}

export interface MolBlockResponse {
  mol_block: string;
}

export interface HitMolBlock {
  name: string;
  mol_block: string;
}
