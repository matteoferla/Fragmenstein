"use client";

import { useState, useEffect } from "react";
import { useParams } from "next/navigation";
import { Button } from "primereact/button";
import { InputNumber } from "primereact/inputnumber";
import { InputText } from "primereact/inputtext";
import { Checkbox } from "primereact/checkbox";
import { TabView, TabPanel } from "primereact/tabview";
import { MolViewer3D } from "@/components/viewer/MolViewer3D";
import { useSessionStore } from "@/stores/sessionStore";
import * as api from "@/services/api";
import type { MonsterResult } from "@/services/types";

export default function MonsterPage() {
  const params = useParams();
  const sessionId = params.id as string;
  const { hits } = useSessionStore();

  const [selectedHits, setSelectedHits] = useState<Set<string>>(new Set());
  const [joiningCutoff, setJoiningCutoff] = useState(5.0);
  const [placeSmiles, setPlaceSmiles] = useState("");
  const [running, setRunning] = useState(false);
  const [result, setResult] = useState<MonsterResult | null>(null);

  useEffect(() => {
    if (hits.length > 0 && selectedHits.size === 0) {
      setSelectedHits(new Set(hits.map(h => h.name)));
    }
  }, [hits, selectedHits.size]);

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
      const res = await api.monsterCombine(sessionId, {
        joining_cutoff: joiningCutoff,
        hit_names: Array.from(selectedHits),
      });
      setResult(res);
    } catch (e: unknown) {
      setResult({ name: "error", smiles: null, num_atoms: null, mol_block: null, error: e instanceof Error ? e.message : "Failed" });
    }
    setRunning(false);
  };

  const handlePlace = async () => {
    if (!placeSmiles) return;
    setRunning(true); setResult(null);
    try {
      const res = await api.monsterPlace(sessionId, {
        smiles: placeSmiles,
        hit_names: Array.from(selectedHits),
      });
      setResult(res);
    } catch (e: unknown) {
      setResult({ name: "error", smiles: null, num_atoms: null, mol_block: null, error: e instanceof Error ? e.message : "Failed" });
    }
    setRunning(false);
  };

  return (
    <div className="max-w-4xl">
      <div className="flex items-center gap-3 mb-6">
        <div className="w-8 h-8 rounded-lg flex items-center justify-center bg-amber-50 border border-amber-200 text-amber-600">
          <i className="pi pi-bolt text-sm" />
        </div>
        <div>
          <h2 className="text-lg font-bold text-slate-800">Monster (Chemistry Only)</h2>
          <p className="text-xs text-slate-400">Combine or place molecules without protein context — fast prototyping</p>
        </div>
      </div>

      {/* Hit selection */}
      {hits.length > 0 && (
        <div className="panel p-5 mb-4">
          <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500 mb-3">
            Select Hits ({selectedHits.size}/{hits.length})
          </h3>
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

      <TabView>
        <TabPanel header="Combine">
          <div className="panel p-5 mt-4">
            <div className="flex items-center gap-4 mb-4">
              <div className="flex flex-col gap-1">
                <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">Joining Cutoff (A)</label>
                <InputNumber value={joiningCutoff} onValueChange={(e) => setJoiningCutoff(e.value ?? 5.0)} min={1} max={20} step={0.5} />
              </div>
            </div>
            <Button label="Combine" icon="pi pi-play" size="small" onClick={handleCombine} loading={running} disabled={selectedHits.size < 2} />
          </div>
        </TabPanel>
        <TabPanel header="Place">
          <div className="panel p-5 mt-4">
            <div className="flex flex-col gap-1 mb-4">
              <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">SMILES to Place</label>
              <InputText value={placeSmiles} onChange={(e) => setPlaceSmiles(e.target.value)} placeholder="e.g. c1ccccc1" />
            </div>
            <Button label="Place" icon="pi pi-play" size="small" onClick={handlePlace} loading={running} disabled={!placeSmiles || selectedHits.size === 0} />
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
              {result.mol_block && (
                <MolViewer3D
                  models={[{ data: result.mol_block, format: "mol" }]}
                  height="350px"
                />
              )}
              <div className="grid grid-cols-3 gap-2">
                <div className="stat-card">
                  <div className="stat-label">Name</div>
                  <div className="stat-value">{result.name}</div>
                </div>
                <div className="stat-card">
                  <div className="stat-label">Heavy Atoms</div>
                  <div className="stat-value text-teal-700">{result.num_atoms ?? "-"}</div>
                </div>
                <div className="stat-card col-span-1">
                  <div className="stat-label">SMILES</div>
                  <div className="text-[10px] font-mono break-all mt-1 text-slate-500">{result.smiles || "-"}</div>
                </div>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
