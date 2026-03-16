"use client";

import { useState, useCallback, useEffect } from "react";
import { useParams, useRouter } from "next/navigation";
import { Button } from "primereact/button";
import { InputNumber } from "primereact/inputnumber";
import { Dropdown } from "primereact/dropdown";
import { JobProgress } from "@/components/jobs/JobProgress";
import { SimilarsTable, type SimilarRow } from "@/components/results/SimilarsTable";
import { DownloadPanel } from "@/components/results/DownloadPanel";
import { useSessionStore } from "@/stores/sessionStore";
import { API_BASE_URL } from "@/lib/constants";
import * as api from "@/services/api";

export default function SimilarsPage() {
  const params = useParams();
  const router = useRouter();
  const sessionId = params.id as string;
  const { combineJobId, similarsJobId, setSimilarsJobId } = useSessionStore();

  const [config, setConfig] = useState({ top_n: 100, dist: 25, length: 200, db: "REAL_dataset", outcome_filter: "acceptable" });
  const [running, setRunning] = useState(false);
  const [results, setResults] = useState<SimilarRow[]>([]);
  const [selectedRow, setSelectedRow] = useState<SimilarRow | null>(null);

  const handleStart = async () => {
    const currentCombineJobId = useSessionStore.getState().combineJobId;
    if (!currentCombineJobId) return;
    setRunning(true); setResults([]);
    try {
      const { job_id } = await api.startSimilars(sessionId, {
        combine_job_id: currentCombineJobId,
        ...config,
      });
      setSimilarsJobId(job_id);
    } catch {
      setRunning(false);
    }
  };

  const handleComplete = useCallback(async () => {
    if (!similarsJobId) return;
    try {
      const res = await api.getJobResults(similarsJobId);
      setResults(res.results as unknown as SimilarRow[]);
    } catch {}
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
    <div className="max-w-7xl">
      <div className="flex items-center gap-3 mb-6">
        <div className="w-8 h-8 rounded-lg flex items-center justify-center bg-indigo-50 border border-indigo-200 text-indigo-500">
          <i className="pi pi-search text-sm" />
        </div>
        <div>
          <h2 className="text-lg font-bold text-slate-800">Find Purchasable Analogs</h2>
          <p className="text-xs text-slate-400">Search SmallWorld for commercially available compounds similar to your mergers</p>
        </div>
      </div>

      {/* Config */}
      {!running && results.length === 0 && (
        <div className="panel p-6 mb-6">
          <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500 mb-4">SmallWorld Search Parameters</h3>
          <div className="grid grid-cols-2 gap-4">
            <div className="flex flex-col gap-1.5">
              <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Top N Mergers to Query</label>
              <InputNumber value={config.top_n} onValueChange={(e) => setConfig(c => ({ ...c, top_n: e.value ?? 100 }))} min={1} max={1000} />
            </div>
            <div className="flex flex-col gap-1.5">
              <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Max Graph Edit Distance</label>
              <InputNumber value={config.dist} onValueChange={(e) => setConfig(c => ({ ...c, dist: e.value ?? 25 }))} min={1} max={100} />
            </div>
            <div className="flex flex-col gap-1.5">
              <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Max Results per Query</label>
              <InputNumber value={config.length} onValueChange={(e) => setConfig(c => ({ ...c, length: e.value ?? 200 }))} min={1} max={1000} />
            </div>
            <div className="flex flex-col gap-1.5">
              <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Combine Outcome Filter</label>
              <Dropdown value={config.outcome_filter} options={["acceptable", "deviant", "equally sized", ""].map(v => ({ label: v || "All", value: v }))} onChange={(e) => setConfig(c => ({ ...c, outcome_filter: e.value }))} />
            </div>
          </div>
          <div className="mt-6">
            <Button label="Search SmallWorld" icon="pi pi-search" size="small" onClick={handleStart} loading={running} />
          </div>
        </div>
      )}

      {/* Progress */}
      {similarsJobId && running && (
        <div className="mb-6"><JobProgress jobId={similarsJobId} onComplete={handleComplete} /></div>
      )}

      {/* Results */}
      {results.length > 0 && (
        <>
          <div className="mb-4 flex items-center justify-between">
            <div className="flex items-center gap-2">
              <span className="font-mono text-lg font-bold text-indigo-600">{results.length}</span>
              <span className="text-xs uppercase tracking-wider text-slate-400">purchasable analogs found</span>
            </div>
            <div className="flex gap-2">
              {similarsJobId && <DownloadPanel jobId={similarsJobId} />}
              <Button label="New Search" icon="pi pi-refresh" severity="secondary" size="small" onClick={() => { setResults([]); setSelectedRow(null); }} />
            </div>
          </div>

          <div className="grid grid-cols-3 gap-5">
            {/* Table */}
            <div className="col-span-2">
              <SimilarsTable results={results} onRowSelect={setSelectedRow} selectedRow={selectedRow} />
            </div>

            {/* Detail Panel */}
            <div>
              {selectedRow ? (
                <div className="panel p-4 space-y-4 sticky top-4">
                  <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500">Analog Detail</h3>

                  {/* 2D Structure */}
                  {selectedRow.smiles && (
                    <div className="flex justify-center bg-white rounded-lg border border-slate-200 p-2">
                      <img
                        src={`${API_BASE_URL}/api/depict?smiles=${encodeURIComponent(selectedRow.smiles)}&width=300&height=220`}
                        alt={selectedRow.smiles}
                        style={{ maxWidth: "100%", height: "auto" }}
                      />
                    </div>
                  )}

                  {/* Metrics */}
                  <div className="grid grid-cols-2 gap-2">
                    <div className="stat-card">
                      <div className="stat-label">Vendor ID</div>
                      <div className="stat-value text-sm">{selectedRow.name || "-"}</div>
                    </div>
                    <div className="stat-card">
                      <div className="stat-label">Topo Distance</div>
                      <div className="stat-value text-teal-700">{selectedRow.topodist ?? "-"}</div>
                    </div>
                    <div className="stat-card">
                      <div className="stat-label">ECFP4 Similarity</div>
                      <div className="stat-value">{selectedRow.ecfp4?.toFixed(3) ?? "-"}</div>
                    </div>
                    <div className="stat-card">
                      <div className="stat-label">Tanimoto</div>
                      <div className="stat-value">{selectedRow.daylight?.toFixed(3) ?? "-"}</div>
                    </div>
                  </div>

                  {/* SMILES */}
                  <div className="stat-card">
                    <div className="stat-label">SMILES</div>
                    <div className="text-[10px] font-mono break-all mt-1 text-slate-500">
                      {selectedRow.smiles || "-"}
                    </div>
                  </div>

                  {/* Query merger */}
                  {selectedRow.query_smiles && (
                    <div className="stat-card">
                      <div className="stat-label">Source Merger</div>
                      <div className="flex items-center gap-2 mt-1">
                        <img
                          src={`${API_BASE_URL}/api/depict?smiles=${encodeURIComponent(selectedRow.query_smiles)}&width=150&height=100`}
                          alt="query"
                          className="rounded border border-slate-100"
                          style={{ width: 75, height: 50, objectFit: "contain", background: "#fff" }}
                        />
                        <span className="text-[9px] font-mono text-slate-400 break-all">
                          {selectedRow.query_smiles}
                        </span>
                      </div>
                    </div>
                  )}
                </div>
              ) : (
                <div className="panel p-8 flex flex-col items-center justify-center text-center" style={{ minHeight: "300px" }}>
                  <i className="pi pi-eye text-2xl mb-3 text-slate-300" />
                  <p className="text-xs text-slate-400">Select an analog to view its structure and details</p>
                </div>
              )}
            </div>
          </div>

          <div className="mt-6 flex justify-end">
            <Button label="Place Analogs" icon="pi pi-arrow-right" iconPos="right" size="small" onClick={() => router.push(`/sessions/${sessionId}/place`)} />
          </div>
        </>
      )}
    </div>
  );
}
