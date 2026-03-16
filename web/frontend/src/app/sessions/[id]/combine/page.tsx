"use client";

import { useState, useCallback, useEffect } from "react";
import { useParams, useRouter } from "next/navigation";
import { Button } from "primereact/button";
import { ConfigForm } from "@/components/common/ConfigForm";
import { JobProgress } from "@/components/jobs/JobProgress";
import { ResultsTable } from "@/components/results/ResultsTable";
import { OutcomeChart } from "@/components/results/OutcomeChart";
import { ResultDetail } from "@/components/results/ResultDetail";
import { DownloadPanel } from "@/components/results/DownloadPanel";
import { useSessionStore } from "@/stores/sessionStore";
import * as api from "@/services/api";
import type { CombineRequest, ResultRow } from "@/services/types";

const COMBINE_FIELDS = [
  { key: "victor_type", label: "Victor Type", type: "select" as const, options: ["Wictor", "Victor", "Quicktor"] },
  { key: "combination_size", label: "Combination Size", type: "number" as const, min: 2, max: 5 },
  { key: "n_cores", label: "Cores (-1 = all)", type: "number" as const, min: -1, max: 64 },
  { key: "timeout", label: "Timeout (s)", type: "number" as const, min: 30, max: 3600 },
  { key: "joining_cutoff", label: "Joining Cutoff (A)", type: "number" as const, min: 1, max: 20, step: 0.5 },
  { key: "permute", label: "Permute", type: "checkbox" as const },
  { key: "covalent_resi", label: "Covalent Residue", type: "text" as const },
];

export default function CombinePage() {
  const params = useParams();
  const router = useRouter();
  const sessionId = params.id as string;
  const { combineJobId, setCombineJobId } = useSessionStore();

  const [config, setConfig] = useState<CombineRequest>({
    victor_type: "Wictor", combination_size: 2, n_cores: -1,
    timeout: 240, joining_cutoff: 5.0, permute: true, covalent_resi: null,
  });
  const [running, setRunning] = useState(false);
  const [results, setResults] = useState<ResultRow[]>([]);
  const [selectedRow, setSelectedRow] = useState<ResultRow | null>(null);

  const handleStart = async () => {
    setRunning(true); setResults([]);
    try { const { job_id } = await api.startCombine(sessionId, config); setCombineJobId(job_id); }
    catch { setRunning(false); }
  };

  const handleComplete = useCallback(async () => {
    if (!combineJobId) return;
    try { const res = await api.getJobResults(combineJobId); setResults(res.results); } catch {}
    setRunning(false);
  }, [combineJobId]);

  useEffect(() => {
    if (combineJobId && results.length === 0) {
      api.getJobStatus(combineJobId).then((job) => {
        if (job.status === "completed") handleComplete();
        else if (job.status === "running") setRunning(true);
      }).catch(() => {});
    }
  }, [combineJobId, results.length, handleComplete]);

  return (
    <div className="max-w-6xl">
      <div className="flex items-center gap-3 mb-6">
        <div className="w-8 h-8 rounded-lg flex items-center justify-center bg-teal-50 border border-teal-200 text-teal-600">
          <i className="pi pi-sitemap text-sm" />
        </div>
        <div>
          <h2 className="text-lg font-bold text-slate-800">Combine Fragments</h2>
          <p className="text-xs text-slate-400">Merge hit fragments into novel ligand candidates</p>
        </div>
      </div>

      {!running && results.length === 0 && (
        <div className="panel p-6 mb-6">
          <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500 mb-4">Parameters</h3>
          <ConfigForm fields={COMBINE_FIELDS} values={config as unknown as Record<string, unknown>} onChange={(k, v) => setConfig(prev => ({ ...prev, [k]: v }))} />
          <div className="mt-6">
            <Button label="Run Combine" icon="pi pi-play" size="small" onClick={handleStart} loading={running} />
          </div>
        </div>
      )}

      {combineJobId && running && (
        <div className="mb-6"><JobProgress jobId={combineJobId} onComplete={handleComplete} /></div>
      )}

      {results.length > 0 && (
        <>
          <div className="mb-4 flex items-center justify-between">
            <div className="flex items-center gap-2">
              <span className="font-mono text-lg font-bold text-teal-600">{results.length}</span>
              <span className="text-xs uppercase tracking-wider text-slate-400">results</span>
            </div>
            <div className="flex gap-2">
              <DownloadPanel jobId={combineJobId!} />
              <Button label="New Run" icon="pi pi-refresh" severity="secondary" size="small" onClick={() => { setResults([]); setSelectedRow(null); }} />
            </div>
          </div>

          <div className="mb-5"><OutcomeChart results={results} /></div>

          <div className="grid grid-cols-3 gap-5">
            <div className="col-span-2">
              <ResultsTable results={results} onRowSelect={setSelectedRow} selectedRow={selectedRow} />
            </div>
            <div>
              {selectedRow ? (
                <ResultDetail result={selectedRow} jobId={combineJobId!} sessionId={sessionId} />
              ) : (
                <div className="panel p-8 flex flex-col items-center justify-center text-center" style={{ minHeight: "300px" }}>
                  <i className="pi pi-eye text-2xl mb-3 text-slate-300" />
                  <p className="text-xs text-slate-400">Select a result row to inspect 3D structure</p>
                </div>
              )}
            </div>
          </div>

          <div className="mt-6 flex justify-end">
            <Button label="Find Similars" icon="pi pi-arrow-right" iconPos="right" size="small" onClick={() => router.push(`/sessions/${sessionId}/similars`)} />
          </div>
        </>
      )}
    </div>
  );
}
