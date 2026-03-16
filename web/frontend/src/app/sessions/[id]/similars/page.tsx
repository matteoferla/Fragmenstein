"use client";

import { useState, useCallback, useEffect } from "react";
import { useParams, useRouter } from "next/navigation";
import { Button } from "primereact/button";
import { InputNumber } from "primereact/inputnumber";
import { Dropdown } from "primereact/dropdown";
import { JobProgress } from "@/components/jobs/JobProgress";
import { ResultsTable } from "@/components/results/ResultsTable";
import { useSessionStore } from "@/stores/sessionStore";
import * as api from "@/services/api";
import type { ResultRow } from "@/services/types";

export default function SimilarsPage() {
  const params = useParams();
  const router = useRouter();
  const sessionId = params.id as string;
  const { combineJobId, similarsJobId, setSimilarsJobId } = useSessionStore();

  const [config, setConfig] = useState({ top_n: 100, dist: 25, length: 200, db: "REAL_dataset", outcome_filter: "acceptable" });
  const [running, setRunning] = useState(false);
  const [results, setResults] = useState<ResultRow[]>([]);

  const handleStart = async () => {
    if (!combineJobId) return;
    setRunning(true); setResults([]);
    try { const { job_id } = await api.startSimilars(sessionId, { combine_job_id: combineJobId, ...config }); setSimilarsJobId(job_id); }
    catch { setRunning(false); }
  };

  const handleComplete = useCallback(async () => {
    if (!similarsJobId) return;
    try { const res = await api.getJobResults(similarsJobId); setResults(res.results); } catch {}
    setRunning(false);
  }, [similarsJobId]);

  useEffect(() => {
    if (similarsJobId && results.length === 0) {
      api.getJobStatus(similarsJobId).then((job) => {
        if (job.status === "completed") handleComplete();
        else if (job.status === "running") setRunning(true);
      }).catch(() => {});
    }
  }, [similarsJobId, results.length, handleComplete]);

  if (!combineJobId) {
    return (
      <div className="max-w-4xl">
        <div className="panel p-8 text-center">
          <i className="pi pi-info-circle text-xl mb-3 text-amber-500" />
          <h2 className="text-lg font-bold mb-2 text-slate-800">Combine Step Required</h2>
          <p className="text-sm mb-4 text-slate-500">Complete the Combine step first to search for purchasable analogs.</p>
          <Button label="Go to Combine" icon="pi pi-arrow-left" size="small" onClick={() => router.push(`/sessions/${sessionId}/combine`)} />
        </div>
      </div>
    );
  }

  return (
    <div className="max-w-6xl">
      <div className="flex items-center gap-3 mb-6">
        <div className="w-8 h-8 rounded-lg flex items-center justify-center bg-indigo-50 border border-indigo-200 text-indigo-500">
          <i className="pi pi-search text-sm" />
        </div>
        <div>
          <h2 className="text-lg font-bold text-slate-800">Find Purchasable Analogs</h2>
          <p className="text-xs text-slate-400">Search SmallWorld for commercially available similar compounds</p>
        </div>
      </div>

      {!running && results.length === 0 && (
        <div className="panel p-6 mb-6">
          <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500 mb-4">SmallWorld Search Parameters</h3>
          <div className="grid grid-cols-2 gap-4">
            <div className="flex flex-col gap-1.5">
              <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Top N Mergers</label>
              <InputNumber value={config.top_n} onValueChange={(e) => setConfig(c => ({ ...c, top_n: e.value ?? 100 }))} min={1} max={1000} />
            </div>
            <div className="flex flex-col gap-1.5">
              <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Distance</label>
              <InputNumber value={config.dist} onValueChange={(e) => setConfig(c => ({ ...c, dist: e.value ?? 25 }))} min={1} max={100} />
            </div>
            <div className="flex flex-col gap-1.5">
              <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Max Results per Query</label>
              <InputNumber value={config.length} onValueChange={(e) => setConfig(c => ({ ...c, length: e.value ?? 200 }))} min={1} max={1000} />
            </div>
            <div className="flex flex-col gap-1.5">
              <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Outcome Filter</label>
              <Dropdown value={config.outcome_filter} options={["acceptable", "deviant", "equally sized", ""].map(v => ({ label: v || "All", value: v }))} onChange={(e) => setConfig(c => ({ ...c, outcome_filter: e.value }))} />
            </div>
          </div>
          <div className="mt-6">
            <Button label="Search SmallWorld" icon="pi pi-search" size="small" onClick={handleStart} loading={running} />
          </div>
        </div>
      )}

      {similarsJobId && running && (
        <div className="mb-6"><JobProgress jobId={similarsJobId} onComplete={handleComplete} /></div>
      )}

      {results.length > 0 && (
        <>
          <div className="flex items-center gap-2 mb-4">
            <span className="font-mono text-lg font-bold text-indigo-600">{results.length}</span>
            <span className="text-xs uppercase tracking-wider text-slate-400">analogs found</span>
          </div>
          <ResultsTable results={results} />
          <div className="mt-6 flex justify-end">
            <Button label="Place Analogs" icon="pi pi-arrow-right" iconPos="right" size="small" onClick={() => router.push(`/sessions/${sessionId}/place`)} />
          </div>
        </>
      )}
    </div>
  );
}
