"use client";

import { useState, useEffect } from "react";
import { useParams } from "next/navigation";
import { Button } from "primereact/button";
import { InputText } from "primereact/inputtext";
import { Dropdown } from "primereact/dropdown";
import { Checkbox } from "primereact/checkbox";
import { TabView, TabPanel } from "primereact/tabview";
import { OverlayViewer } from "@/components/viewer/OverlayViewer";
import { useSessionStore } from "@/stores/sessionStore";
import { VICTOR_TYPES } from "@/lib/constants";
import * as api from "@/services/api";
import type { SingleResult } from "@/services/types";

export default function SingleVictorPage() {
  const params = useParams();
  const sessionId = params.id as string;
  const { hits } = useSessionStore();

  const [selectedHits, setSelectedHits] = useState<Set<string>>(new Set());
  const [victorType, setVictorType] = useState("Wictor");
  const [covalentResi, setCovalentResi] = useState("");
  const [placeSmiles, setPlaceSmiles] = useState("");
  const [running, setRunning] = useState(false);
  const [result, setResult] = useState<SingleResult | null>(null);
  const [proteinPdb, setProteinPdb] = useState<string | undefined>();
  const [hitMolBlocks, setHitMolBlocks] = useState<Array<{ name: string; data: string }>>([]);

  useEffect(() => {
    if (hits.length > 0 && selectedHits.size === 0) {
      setSelectedHits(new Set(hits.slice(0, 2).map(h => h.name)));
    }
  }, [hits, selectedHits.size]);

  useEffect(() => {
    api.getTemplatePdb(sessionId).then(r => setProteinPdb(r.pdb)).catch(() => {});
  }, [sessionId]);

  // Load hit mol blocks for selected hits
  useEffect(() => {
    const names = Array.from(selectedHits);
    if (names.length === 0) { setHitMolBlocks([]); return; }
    Promise.all(names.map(n => api.getHitMolBlock(sessionId, n).then(r => ({ name: n, data: r.mol_block })).catch(() => null)))
      .then(results => setHitMolBlocks(results.filter(Boolean) as Array<{ name: string; data: string }>));
  }, [selectedHits, sessionId]);

  const toggleHit = (name: string) => {
    setSelectedHits(prev => {
      const next = new Set(prev);
      if (next.has(name)) next.delete(name); else next.add(name);
      return next;
    });
  };

  const handleCombine = async () => {
    setRunning(true); setResult(null);
    try {
      const res = await api.singleCombine(sessionId, {
        hit_names: Array.from(selectedHits),
        victor_type: victorType,
        covalent_resi: covalentResi || null,
      });
      setResult(res);
    } catch (e: unknown) {
      setResult({ name: "error", smiles: null, error: e instanceof Error ? e.message : "Failed", mode: null, ddG: null, dG_bound: null, dG_unbound: null, comRMSD: null, N_constrained_atoms: null, N_unconstrained_atoms: null, runtime: null, mol_block: null, unmin_mol_block: null });
    }
    setRunning(false);
  };

  const handlePlace = async () => {
    if (!placeSmiles) return;
    setRunning(true); setResult(null);
    try {
      const res = await api.singlePlace(sessionId, {
        smiles: placeSmiles,
        hit_names: Array.from(selectedHits),
        victor_type: victorType,
        covalent_resi: covalentResi || null,
      });
      setResult(res);
    } catch (e: unknown) {
      setResult({ name: "error", smiles: null, error: e instanceof Error ? e.message : "Failed", mode: null, ddG: null, dG_bound: null, dG_unbound: null, comRMSD: null, N_constrained_atoms: null, N_unconstrained_atoms: null, runtime: null, mol_block: null, unmin_mol_block: null });
    }
    setRunning(false);
  };

  return (
    <div className="max-w-4xl">
      <div className="flex items-center gap-3 mb-6">
        <div className="w-8 h-8 rounded-lg flex items-center justify-center bg-violet-50 border border-violet-200 text-violet-600">
          <i className="pi pi-wrench text-sm" />
        </div>
        <div>
          <h2 className="text-lg font-bold text-slate-800">Single Victor</h2>
          <p className="text-xs text-slate-400">Run one combine or place with full energy scoring — for debugging and exploration</p>
        </div>
      </div>

      {/* Hit selection */}
      {hits.length > 0 && (
        <div className="panel p-5 mb-4">
          <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500 mb-3">Select Hits</h3>
          <div className="flex flex-wrap gap-2">
            {hits.map(h => (
              <label key={h.name} className={`flex items-center gap-1.5 px-2.5 py-1.5 rounded-lg border cursor-pointer text-xs transition-colors ${
                selectedHits.has(h.name) ? "bg-teal-50 border-teal-200 text-teal-700" : "bg-slate-50 border-slate-200 text-slate-400"
              }`}>
                <Checkbox checked={selectedHits.has(h.name)} onChange={() => toggleHit(h.name)} />
                <span className="font-mono">{h.name}</span>
              </label>
            ))}
          </div>
        </div>
      )}

      {/* Config */}
      <div className="panel p-5 mb-4">
        <div className="grid grid-cols-2 gap-4">
          <div className="flex flex-col gap-1.5">
            <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Victor Type</label>
            <Dropdown value={victorType} options={VICTOR_TYPES.map(v => ({ label: v, value: v }))} onChange={e => setVictorType(e.value)} />
          </div>
          <div className="flex flex-col gap-1.5">
            <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Covalent Residue</label>
            <InputText value={covalentResi} onChange={e => setCovalentResi(e.target.value)} placeholder="e.g. 145A" />
          </div>
        </div>
      </div>

      <TabView>
        <TabPanel header="Combine">
          <div className="mt-4">
            <Button label="Run Single Combine" icon="pi pi-play" size="small" onClick={handleCombine} loading={running} disabled={selectedHits.size < 2} />
          </div>
        </TabPanel>
        <TabPanel header="Place">
          <div className="mt-4 space-y-3">
            <div className="flex flex-col gap-1.5">
              <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">SMILES to Place</label>
              <InputText value={placeSmiles} onChange={e => setPlaceSmiles(e.target.value)} placeholder="e.g. c1ccc(CC(=O)N)cc1" />
            </div>
            <Button label="Run Single Place" icon="pi pi-play" size="small" onClick={handlePlace} loading={running} disabled={!placeSmiles || selectedHits.size === 0} />
          </div>
        </TabPanel>
      </TabView>

      {/* Result */}
      {result && (
        <div className="panel p-5 mt-6">
          <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500 mb-3">Result</h3>
          {result.error ? (
            <p className="text-sm text-red-600">{result.error}</p>
          ) : (
            <div className="space-y-4">
              <OverlayViewer
                proteinPdb={proteinPdb}
                hitMolBlocks={hitMolBlocks}
                resultMolBlock={result.mol_block ?? undefined}
                height="400px"
              />
              <div className="grid grid-cols-3 gap-2">
                <div className="stat-card">
                  <div className="stat-label">ddG</div>
                  <div className="stat-value text-teal-700">{result.ddG?.toFixed(2) ?? "-"} <span className="text-[10px] text-slate-400">kcal/mol</span></div>
                </div>
                <div className="stat-card">
                  <div className="stat-label">comRMSD</div>
                  <div className="stat-value">{result.comRMSD?.toFixed(2) ?? "-"} <span className="text-[10px] text-slate-400">A</span></div>
                </div>
                <div className="stat-card">
                  <div className="stat-label">Runtime</div>
                  <div className="stat-value">{result.runtime?.toFixed(1) ?? "-"} <span className="text-[10px] text-slate-400">s</span></div>
                </div>
              </div>
              <div className="stat-card">
                <div className="stat-label">SMILES</div>
                <div className="text-[10px] font-mono break-all mt-1 text-slate-500">{result.smiles || "-"}</div>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
