/** API client for the Fragmenstein backend. */

import { API_BASE_URL } from "@/lib/constants";
import type {
  CombineRequest,
  HitsResponse,
  JobStatus,
  MolBlockResponse,
  PlaceRequest,
  ResultsResponse,
  SessionResponse,
  SimilarsRequest,
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

export async function uploadHits(sessionId: string, files: File[]): Promise<{ filename: string; message: string }> {
  const formData = new FormData();
  files.forEach((f) => formData.append("files", f));
  const res = await fetch(`${API_BASE_URL}/api/sessions/${sessionId}/hits`, {
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

export async function getJobResults(jobId: string): Promise<ResultsResponse> {
  return request(`/api/jobs/${jobId}/results`);
}

export async function getResultMol(jobId: string, idx: number, molType: string = "minimized"): Promise<MolBlockResponse> {
  return request(`/api/jobs/${jobId}/results/${idx}/mol?mol_type=${molType}`);
}

export function getResultsDownloadUrl(jobId: string, format: string = "csv"): string {
  return `${API_BASE_URL}/api/jobs/${jobId}/results/download?format=${format}`;
}

// Molecules
export async function getHitMolBlock(sessionId: string, hitName: string): Promise<MolBlockResponse> {
  return request(`/api/sessions/${sessionId}/hits/${hitName}/mol`);
}

export async function getAllHitMolBlocks(sessionId: string): Promise<{ hits: HitMolBlock[] }> {
  return request(`/api/sessions/${sessionId}/hits/molblocks`);
}
