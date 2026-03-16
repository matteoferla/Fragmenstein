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
import type { PlaceRequest, ResultRow } from "@/services/types";

const PLACE_FIELDS = [
  { key: "victor_type", label: "Victor Type", type: "select" as const, options: ["Wictor", "Victor", "Quicktor"] },
  { key: "n_cores", label: "Cores (-1 = all)", type: "number" as const, min: -1, max: 64 },
  { key: "timeout", label: "Timeout (s)", type: "number" as const, min: 30, max: 3600 },
  { key: "covalent_resi", label: "Covalent Residue", type: "text" as const },
];

export default function PlacePage() {
  const params = useParams();
  const router = useRouter();
  const sessionId = params.id as string;
  const { similarsJobId, placeJobId, setPlaceJobId } = useSessionStore();

  const [config, setConfig] = useState<PlaceRequest>({ victor_type: "Wictor", n_cores: -1, timeout: 240, covalent_resi: null, source_job_id: similarsJobId });
  const [running, setRunning] = useState(false);
  const [results, setResults] = useState<ResultRow[]>([]);
  const [selectedRow, setSelectedRow] = useState<ResultRow | null>(null);

  useEffect(() => { if (similarsJobId) setConfig(c => ({ ...c, source_job_id: similarsJobId })); }, [similarsJobId]);

  const handleStart = async () => {
    if (!config.source_job_id) return;
    setRunning(true); setResults([]);
    try { const { job_id } = await api.startPlace(sessionId, config); setPlaceJobId(job_id); }
    catch { setRunning(false); }
  };

  const handleComplete = useCallback(async () => {
    if (!placeJobId) return;
    try { const res = await api.getJobResults(placeJobId); setResults(res.results); } catch {}
    setRunning(false);
  }, [placeJobId]);

  useEffect(() => {
    if (placeJobId && results.length === 0) {
      api.getJobStatus(placeJobId).then((job) => {
        if (job.status === "completed") handleComplete();
        else if (job.status === "running") setRunning(true);
      }).catch(() => {});
    }
  }, [placeJobId, results.length, handleComplete]);

  if (!similarsJobId) {
    return (
      <div className="max-w-4xl">
        <div className="panel p-8 text-center">
          <i className="pi pi-info-circle text-xl mb-3 text-amber-500" />
          <h2 className="text-lg font-bold mb-2 text-slate-800">Similars Step Required</h2>
          <p className="text-sm mb-4 text-slate-500">Complete the Similars step first to place analog molecules.</p>
          <Button label="Go to Similars" icon="pi pi-arrow-left" size="small" onClick={() => router.push(`/sessions/${sessionId}/similars`)} />
        </div>
      </div>
    );
  }

  return (
    <div className="max-w-6xl">
      <div className="flex items-center gap-3 mb-6">
        <div className="w-8 h-8 rounded-lg flex items-center justify-center bg-emerald-50 border border-emerald-200 text-emerald-600">
          <i className="pi pi-map-marker text-sm" />
        </div>
        <div>
          <h2 className="text-lg font-bold text-slate-800">Place Analogs</h2>
          <p className="text-xs text-slate-400">Dock analog molecules into the binding site</p>
        </div>
      </div>

      {!running && results.length === 0 && (
        <div className="panel p-6 mb-6">
          <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500 mb-4">Placement Parameters</h3>
          <ConfigForm fields={PLACE_FIELDS} values={config as unknown as Record<string, unknown>} onChange={(k, v) => setConfig(prev => ({ ...prev, [k]: v }))} />
          <div className="mt-6">
            <Button label="Run Placement" icon="pi pi-play" size="small" onClick={handleStart} loading={running} />
          </div>
        </div>
      )}

      {placeJobId && running && (
        <div className="mb-6"><JobProgress jobId={placeJobId} onComplete={handleComplete} /></div>
      )}

      {results.length > 0 && (
        <>
          <div className="mb-4 flex items-center justify-between">
            <div className="flex items-center gap-2">
              <span className="font-mono text-lg font-bold text-emerald-600">{results.length}</span>
              <span className="text-xs uppercase tracking-wider text-slate-400">placements</span>
            </div>
            <DownloadPanel jobId={placeJobId!} />
          </div>
          <div className="mb-5"><OutcomeChart results={results} /></div>
          <div className="grid grid-cols-3 gap-5">
            <div className="col-span-2">
              <ResultsTable results={results} onRowSelect={setSelectedRow} selectedRow={selectedRow} />
            </div>
            <div>
              {selectedRow ? (
                <ResultDetail result={selectedRow} jobId={placeJobId!} sessionId={sessionId} />
              ) : (
                <div className="panel p-8 flex flex-col items-center justify-center text-center" style={{ minHeight: "300px" }}>
                  <i className="pi pi-eye text-2xl mb-3 text-slate-300" />
                  <p className="text-xs text-slate-400">Select a result row to inspect 3D structure</p>
                </div>
              )}
            </div>
          </div>
          <div className="mt-6 flex justify-end">
            <Button label="View All Results" icon="pi pi-arrow-right" iconPos="right" size="small" onClick={() => router.push(`/sessions/${sessionId}/results`)} />
          </div>
        </>
      )}
    </div>
  );
}
