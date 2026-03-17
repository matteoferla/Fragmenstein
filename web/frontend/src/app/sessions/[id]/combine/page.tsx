"use client";

import { useState, useCallback, useEffect } from "react";
import { useParams, useRouter } from "next/navigation";
import { Button } from "primereact/button";
import { Checkbox } from "primereact/checkbox";
import { ConfigForm } from "@/components/common/ConfigForm";
import { JobProgress } from "@/components/jobs/JobProgress";
import { ResultsTable } from "@/components/results/ResultsTable";
import { OutcomeChart } from "@/components/results/OutcomeChart";
import { ResultDetail } from "@/components/results/ResultDetail";
import { DownloadPanel } from "@/components/results/DownloadPanel";
import { MetricsLegend } from "@/components/results/MetricsLegend";
import { useSessionStore } from "@/stores/sessionStore";
import { VICTOR_TYPES } from "@/lib/constants";
import * as api from "@/services/api";
import type { CombineRequest, ResultRow } from "@/services/types";

const COMBINE_FIELDS = [
  {
    key: "victor_type", label: "Victor Type", type: "select" as const, options: VICTOR_TYPES,
    optionDescs: {
      Wictor: "RDKit-only minimization. Fast (~20s), no PyRosetta needed.",
      Victor: "Full PyRosetta energy scoring. Slow (~60s) but most accurate.",
      Quicktor: "Quick PyRosetta mode. Medium speed, strict MCS matching.",
      OpenVictor: "OpenMM minimization. GPU-capable, free alternative to PyRosetta.",
    },
  },
  { key: "combination_size", label: "Combination Size", type: "number" as const, min: 2, max: 5, description: "Number of fragments to merge at once. 2 = pairwise, 3 = triplets." },
  { key: "n_cores", label: "CPU Cores", type: "number" as const, min: -1, max: 64, description: "-1 uses all available cores. Positive number = exact core count." },
  { key: "timeout", label: "Timeout (s)", type: "number" as const, min: 30, max: 3600, description: "Max seconds per combination. Longer = more results but slower." },
  { key: "joining_cutoff", label: "Joining Cutoff (A)", type: "number" as const, min: 1, max: 20, step: 0.5, description: "Max distance (Angstroms) to link disconnected fragments with a bond." },
  {
    key: "warhead_harmonisation", label: "Warhead Mode", type: "select" as const,
    options: ["first", "keep", "strip", "acrylamide", "chloroacetamide", "nitrile", "vinylsulfonamide", "bromoalkyne"],
    optionDescs: {
      first: "Use the warhead from the first hit fragment.",
      keep: "Keep all warheads as-is from each fragment.",
      strip: "Remove all warheads — produce non-covalent mergers only.",
      acrylamide: "Force acrylamide warhead on the merger.",
      chloroacetamide: "Force chloroacetamide warhead.",
      nitrile: "Force nitrile warhead.",
      vinylsulfonamide: "Force vinylsulfonamide warhead.",
      bromoalkyne: "Force bromoalkyne warhead.",
    },
  },
  { key: "covalent_resi", label: "Covalent Residue", type: "text" as const, description: "PDB residue for covalent attachment, e.g. '145A'. Leave empty for non-covalent." },
  { key: "permute", label: "Permute", type: "checkbox" as const, description: "Generate all orderings (A+B and B+A). More results but slower." },
  { key: "quick_reanimation", label: "Quick Reanimation", type: "checkbox" as const, description: "Faster but less thorough energy minimization." },
  { key: "run_plip", label: "PLIP Analysis", type: "checkbox" as const, description: "Run protein-ligand interaction analysis (H-bonds, hydrophobic contacts)." },
];

export default function CombinePage() {
  const params = useParams();
  const router = useRouter();
  const sessionId = params.id as string;
  const { combineJobId, setCombineJobId, hits } = useSessionStore();

  const [config, setConfig] = useState<CombineRequest>({
    victor_type: "Wictor", combination_size: 2, n_cores: -1,
    timeout: 240, joining_cutoff: 5.0, permute: true,
    warhead_harmonisation: "first", quick_reanimation: false,
    run_plip: false, covalent_resi: null, hit_names: null,
  });
  const [running, setRunning] = useState(false);
  const [results, setResults] = useState<ResultRow[]>([]);
  const [selectedRow, setSelectedRow] = useState<ResultRow | null>(null);
  const [selectedHits, setSelectedHits] = useState<Set<string>>(new Set());

  // Initialize hit selection from session hits
  useEffect(() => {
    if (hits.length > 0 && selectedHits.size === 0) {
      setSelectedHits(new Set(hits.map(h => h.name)));
    }
  }, [hits, selectedHits.size]);

  const toggleHit = (name: string) => {
    setSelectedHits(prev => {
      const next = new Set(prev);
      if (next.has(name)) next.delete(name);
      else next.add(name);
      return next;
    });
  };

  const handleStart = async () => {
    const hitNames = selectedHits.size === hits.length ? null : Array.from(selectedHits);
    const finalConfig = { ...config, hit_names: hitNames };
    setRunning(true); setResults([]);
    try { const { job_id } = await api.startCombine(sessionId, finalConfig); setCombineJobId(job_id); }
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
        <div className="w-8 h-8 rounded-lg flex items-center justify-center bg-blue-50 border border-blue-200 text-blue-600">
          <i className="pi pi-sitemap text-sm" />
        </div>
        <div>
          <h2 className="text-lg font-bold text-slate-800">Combine Fragments</h2>
          <p className="text-xs text-slate-400">Merge hit fragments into novel ligand candidates</p>
        </div>
      </div>

      {!running && results.length === 0 && (
        <>
          {/* Hit Selection */}
          {hits.length > 0 && (
            <div className="panel p-5 mb-4">
              <div className="flex items-center justify-between mb-3">
                <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500">
                  Select Hits ({selectedHits.size}/{hits.length})
                </h3>
                <div className="flex gap-2">
                  <Button label="All" size="small" severity="secondary" className="!text-[10px] !py-1 !px-2" onClick={() => setSelectedHits(new Set(hits.map(h => h.name)))} />
                  <Button label="None" size="small" severity="secondary" className="!text-[10px] !py-1 !px-2" onClick={() => setSelectedHits(new Set())} />
                </div>
              </div>
              <div className="flex flex-wrap gap-2">
                {hits.map((h, i) => (
                  <label key={`${h.name}-${i}`} className={`flex items-center gap-1.5 px-2.5 py-1.5 rounded-lg border cursor-pointer text-xs transition-colors ${
                    selectedHits.has(h.name) ? "bg-blue-50 border-blue-200 text-blue-700" : "bg-slate-50 border-slate-200 text-slate-400"
                  }`}>
                    <Checkbox checked={selectedHits.has(h.name)} onChange={() => toggleHit(h.name)} />
                    <span className="font-mono">{h.name}</span>
                    {h.num_atoms && <span className="text-[10px] text-slate-400">({h.num_atoms} HA)</span>}
                  </label>
                ))}
              </div>
            </div>
          )}

          {/* Parameters */}
          <div className="panel p-6 mb-6">
            <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500 mb-4">Parameters</h3>
            <ConfigForm fields={COMBINE_FIELDS} values={config as unknown as Record<string, unknown>} onChange={(k, v) => setConfig(prev => ({ ...prev, [k]: v }))} />
            <div className="mt-6">
              <Button label="Run Combine" icon="pi pi-play" size="small" onClick={handleStart} loading={running} disabled={selectedHits.size < 2} />
            </div>
          </div>
        </>
      )}

      {combineJobId && running && (
        <div className="mb-6">
          <JobProgress
            jobId={combineJobId}
            onComplete={handleComplete}
            onCancel={() => setRunning(false)}
            onRerun={() => { setRunning(false); setResults([]); setSelectedRow(null); setCombineJobId(null); }}
          />
        </div>
      )}

      {results.length > 0 && (
        <>
          <div className="mb-4 flex items-center justify-between">
            <div className="flex items-center gap-2">
              <span className="font-mono text-lg font-bold text-blue-600">{results.length}</span>
              <span className="text-xs uppercase tracking-wider text-slate-400">results</span>
            </div>
            <div className="flex gap-2">
              <DownloadPanel jobId={combineJobId!} />
              <Button label="New Run" icon="pi pi-refresh" severity="secondary" size="small" onClick={() => { setResults([]); setSelectedRow(null); setCombineJobId(null); }} />
            </div>
          </div>

          <div className="mb-5"><OutcomeChart results={results} /></div>
          <MetricsLegend />
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
