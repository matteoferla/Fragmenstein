/** API client for the Fragmenstein backend. */

import { API_BASE_URL } from "@/lib/constants";
import type {
  AvailableBackends,
  ChemSpaceRequest,
  CombineRequest,
  HitsResponse,
  JobStatus,
  MolBlockResponse,
  MolPortRequest,
  MonsterCombineRequest,
  MonsterPlaceRequest,
  MonsterResult,
  PlaceRequest,
  ResultsResponse,
  SessionResponse,
  SimilarsRequest,
  SingleCombineRequest,
  TemplatePrepRequest,
  SinglePlaceRequest,
  SingleResult,
  HitMolBlock,
} from "./types";

async function request<T>(path: string, options?: RequestInit): Promise<T> {
  const res = await fetch(`${API_BASE_URL}${path}`, {
    headers: { "Content-Type": "application/json", ...options?.headers },
    ...options,
  });
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`API error ${res.status}: ${body}`);
  }
  return res.json();
}

// System
export async function getSystemInfo(): Promise<Record<string, unknown>> {
  return request("/api/system-info");
}

// Sessions
export async function createSession(name: string = ""): Promise<SessionResponse> {
  return request("/api/sessions", {
    method: "POST",
    body: JSON.stringify({ name }),
  });
}

export async function listSessions(): Promise<SessionResponse[]> {
  return request("/api/sessions");
}

export async function getSession(id: string): Promise<SessionResponse> {
  return request(`/api/sessions/${id}`);
}

export async function deleteSession(id: string): Promise<void> {
  await fetch(`${API_BASE_URL}/api/sessions/${id}`, { method: "DELETE" });
}

export async function getSessionJobs(sessionId: string): Promise<JobStatus[]> {
  return request(`/api/sessions/${sessionId}/jobs`);
}

// Uploads
export async function uploadTemplate(sessionId: string, file: File): Promise<{ filename: string; message: string }> {
  const formData = new FormData();
  formData.append("file", file);
  const res = await fetch(`${API_BASE_URL}/api/sessions/${sessionId}/template`, {
    method: "POST",
    body: formData,
  });
  if (!res.ok) throw new Error(`Upload failed: ${res.status}`);
  return res.json();
}

export async function uploadHits(
  sessionId: string,
  files: File[],
  ligandResn?: string,
  proximityBonding: boolean = true,
): Promise<{ filename: string; message: string }> {
  const formData = new FormData();
  files.forEach((f) => formData.append("files", f));
  const params = new URLSearchParams();
  if (ligandResn) params.set("ligand_resn", ligandResn);
  params.set("proximity_bonding", String(proximityBonding));
  const url = `${API_BASE_URL}/api/sessions/${sessionId}/hits?${params.toString()}`;
  const res = await fetch(url, {
    method: "POST",
    body: formData,
  });
  if (!res.ok) throw new Error(`Upload failed: ${res.status}`);
  return res.json();
}

export async function getHits(sessionId: string): Promise<HitsResponse> {
  return request(`/api/sessions/${sessionId}/hits`);
}

export async function getTemplatePdb(sessionId: string): Promise<{ pdb: string }> {
  return request(`/api/sessions/${sessionId}/template/pdb`);
}

// Template cleaning
export async function cleanTemplate(sessionId: string, removeResidues: string = "HOH"): Promise<{ message: string; removed_count: number }> {
  return request(`/api/sessions/${sessionId}/template/clean?remove_residues=${encodeURIComponent(removeResidues)}`, {
    method: "POST",
  });
}

// Template preparation (PyRosetta)
export async function prepareTemplate(sessionId: string, config: TemplatePrepRequest): Promise<{ job_id: string }> {
  return request(`/api/sessions/${sessionId}/template/prepare`, {
    method: "POST",
    body: JSON.stringify(config),
  });
}

// Demo
export async function loadDemo(sessionId: string, dataset: string, nHits: number = 5): Promise<{ template_filename: string; hit_count: number; message: string }> {
  return request(`/api/sessions/${sessionId}/demo/${dataset}?n_hits=${nHits}`, {
    method: "POST",
  });
}

// Combine
export async function startCombine(sessionId: string, config: CombineRequest): Promise<{ job_id: string }> {
  return request(`/api/sessions/${sessionId}/combine`, {
    method: "POST",
    body: JSON.stringify(config),
  });
}

// Similars
export async function startSimilars(sessionId: string, config: SimilarsRequest): Promise<{ job_id: string }> {
  return request(`/api/sessions/${sessionId}/similars`, {
    method: "POST",
    body: JSON.stringify(config),
  });
}

// Similars — manual paste
export async function manualSmiles(sessionId: string, smilesText: string): Promise<{ job_id: string; valid: number; invalid: number; message: string }> {
  return request(`/api/sessions/${sessionId}/similars/manual`, {
    method: "POST",
    body: JSON.stringify({ smiles_text: smilesText }),
  });
}

// Similars — CSV/Excel upload
export async function uploadSimilars(sessionId: string, file: File): Promise<{ job_id: string; valid: number; invalid: number; message: string }> {
  const formData = new FormData();
  formData.append("file", file);
  const res = await fetch(`${API_BASE_URL}/api/sessions/${sessionId}/similars/upload`, {
    method: "POST",
    body: formData,
  });
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`Upload failed: ${body}`);
  }
  return res.json();
}

// Similars — upload + filter against mergers
export async function uploadAndFilterSimilars(
  sessionId: string,
  file: File,
  topN: number = 200,
  outcomeFilter: string = "acceptable",
): Promise<{ job_id: string; library_size: number; filtered: number; mergers_used: number; invalid: number; message: string }> {
  const formData = new FormData();
  formData.append("file", file);
  const params = new URLSearchParams({ top_n: String(topN), outcome_filter: outcomeFilter });
  const res = await fetch(`${API_BASE_URL}/api/sessions/${sessionId}/similars/upload-filter?${params}`, {
    method: "POST",
    body: formData,
  });
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`Upload failed: ${body}`);
  }
  return res.json();
}

// Similars — PubChem
export async function startPubChem(sessionId: string, config: { combine_job_id?: string | null; top_n?: number; threshold?: number; max_per_query?: number; outcome_filter?: string }): Promise<{ job_id: string }> {
  const { combine_job_id, ...rest } = config;
  const body = combine_job_id ? { combine_job_id, ...rest } : rest;
  return request(`/api/sessions/${sessionId}/similars/pubchem`, {
    method: "POST",
    body: JSON.stringify(body),
  });
}

// Similars — ChemSpace
export async function startChemSpace(sessionId: string, config: ChemSpaceRequest): Promise<{ job_id: string }> {
  const { combine_job_id, ...rest } = config;
  const body = combine_job_id ? { combine_job_id, ...rest } : rest;
  return request(`/api/sessions/${sessionId}/similars/chemspace`, {
    method: "POST",
    body: JSON.stringify(body),
  });
}

// Similars — MolPort
export async function startMolPort(sessionId: string, config: MolPortRequest): Promise<{ job_id: string }> {
  const { combine_job_id, ...rest } = config;
  const body = combine_job_id ? { combine_job_id, ...rest } : rest;
  return request(`/api/sessions/${sessionId}/similars/molport`, {
    method: "POST",
    body: JSON.stringify(body),
  });
}

// Config
export async function getAvailableBackends(): Promise<AvailableBackends> {
  return request("/api/config/available-backends");
}

// Place
export async function startPlace(sessionId: string, config: PlaceRequest): Promise<{ job_id: string }> {
  return request(`/api/sessions/${sessionId}/place`, {
    method: "POST",
    body: JSON.stringify(config),
  });
}

// Jobs
export async function getJobStatus(jobId: string): Promise<JobStatus> {
  return request(`/api/jobs/${jobId}`);
}

export async function cancelJob(jobId: string): Promise<{ status: string; message: string }> {
  return request(`/api/jobs/${jobId}/cancel`, { method: "POST" });
}

export async function getJobResults(jobId: string): Promise<ResultsResponse> {
  return request(`/api/jobs/${jobId}/results`);
}

export async function getResultMol(jobId: string, idx: number, molType: string = "minimized"): Promise<MolBlockResponse> {
  return request(`/api/jobs/${jobId}/results/${idx}/mol?mol_type=${molType}`);
}

export function getResultPdbUrl(jobId: string, idx: number): string {
  return `${API_BASE_URL}/api/jobs/${jobId}/results/${idx}/pdb`;
}

export function getResultsDownloadUrl(jobId: string, format: string = "csv"): string {
  return `${API_BASE_URL}/api/jobs/${jobId}/results/download?format=${format}`;
}

// Monster
export async function monsterCombine(sessionId: string, config: MonsterCombineRequest): Promise<MonsterResult> {
  return request(`/api/sessions/${sessionId}/monster/combine`, { method: "POST", body: JSON.stringify(config) });
}

export async function monsterPlace(sessionId: string, config: MonsterPlaceRequest): Promise<MonsterResult> {
  return request(`/api/sessions/${sessionId}/monster/place`, { method: "POST", body: JSON.stringify(config) });
}

// Single Victor
export async function singleCombine(sessionId: string, config: SingleCombineRequest): Promise<SingleResult> {
  return request(`/api/sessions/${sessionId}/single/combine`, { method: "POST", body: JSON.stringify(config) });
}

export async function singlePlace(sessionId: string, config: SinglePlaceRequest): Promise<SingleResult> {
  return request(`/api/sessions/${sessionId}/single/place`, { method: "POST", body: JSON.stringify(config) });
}

// Molecules
export async function getHitMolBlock(sessionId: string, hitName: string): Promise<MolBlockResponse> {
  return request(`/api/sessions/${sessionId}/hits/${hitName}/mol`);
}

export async function getAllHitMolBlocks(sessionId: string): Promise<{ hits: HitMolBlock[] }> {
  return request(`/api/sessions/${sessionId}/hits/molblocks`);
}
