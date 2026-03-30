"use client";

import { useEffect, useState } from "react";
import { Button } from "primereact/button";
import { OverlayViewer } from "@/components/viewer/OverlayViewer";
import * as api from "@/services/api";
import { getResultPdbUrl } from "@/services/api";
import type { ResultRow } from "@/services/types";
import { OUTCOME_COLORS } from "@/lib/constants";

interface ResultDetailProps {
  result: ResultRow;
  jobId: string;
  sessionId: string;
}

export function ResultDetail({ result, jobId, sessionId }: ResultDetailProps) {
  const [molBlock, setMolBlock] = useState<string | undefined>();
  const [proteinPdb, setProteinPdb] = useState<string | undefined>();
  const [hitMolBlocks, setHitMolBlocks] = useState<Array<{ name: string; data: string }>>([]);

  useEffect(() => {
    api.getResultMol(jobId, result.index).then((res) => setMolBlock(res.mol_block)).catch(() => {});
    api.getTemplatePdb(sessionId).then((res) => setProteinPdb(res.pdb)).catch(() => {});
    if (result.hit_names && result.hit_names.length > 0) {
      Promise.all(
        result.hit_names.map((name) =>
          api.getHitMolBlock(sessionId, name).then((res) => ({ name, data: res.mol_block }))
        )
      ).then(setHitMolBlocks).catch(() => {});
    }
  }, [result, jobId, sessionId]);

  const outcomeColor = OUTCOME_COLORS[result.outcome || ""] || "#94a3b8";

  return (
    <div className="panel p-4 space-y-4">
      <div className="flex items-center justify-between">
        <h3 className="text-sm font-bold font-mono text-slate-800">{result.name}</h3>
        <span
          className="text-[10px] font-bold uppercase tracking-wider font-mono px-2 py-0.5 rounded"
          style={{ background: `${outcomeColor}15`, color: outcomeColor, border: `1px solid ${outcomeColor}30` }}
        >
          {result.outcome || "unknown"}
        </span>
      </div>

      <OverlayViewer proteinPdb={proteinPdb} hitMolBlocks={hitMolBlocks} resultMolBlock={molBlock} height="350px" />

      {/* Download single result PDB */}
      <div className="flex justify-end">
        <a href={getResultPdbUrl(jobId, result.index)} download>
          <Button label="Download PDB" icon="pi pi-download" size="small" severity="secondary" />
        </a>
      </div>

      <div className="grid grid-cols-2 gap-2">
        <div className="stat-card">
          <div className="stat-label">ddG</div>
          <div className="stat-value text-blue-700">
            {result.ddG?.toFixed(2) ?? "-"} <span className="text-[10px] text-slate-400">kcal/mol</span>
          </div>
        </div>
        <div className="stat-card">
          <div className="stat-label">comRMSD</div>
          <div className="stat-value">
            {result.comRMSD?.toFixed(2) ?? "-"} <span className="text-[10px] text-slate-400">A</span>
          </div>
        </div>
        <div className="stat-card">
          <div className="stat-label">Ligand Efficiency</div>
          <div className="stat-value">{result.LE?.toFixed(3) ?? "-"}</div>
        </div>
        <div className="stat-card">
          <div className="stat-label">Heavy Atoms</div>
          <div className="stat-value text-blue-700">{result.N_HA ?? "-"}</div>
        </div>
      </div>

      <div className="stat-card">
        <div className="stat-label">SMILES</div>
        <div className="text-[11px] font-mono break-all mt-1 text-slate-500">{result.smiles || "-"}</div>
      </div>

      {result.hit_names && result.hit_names.length > 0 && (
        <div className="stat-card">
          <div className="stat-label">Source Hits</div>
          <div className="flex gap-1 mt-1 flex-wrap">
            {result.hit_names.map((name, i) => (
              <span key={`${name}-${i}`} className="text-[10px] font-mono px-1.5 py-0.5 rounded bg-emerald-50 text-emerald-700 border border-emerald-200">
                {name}
              </span>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
