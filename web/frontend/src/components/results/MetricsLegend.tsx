"use client";

import { useState } from "react";

const METRICS = [
  { key: "ddG", unit: "kcal/mol", label: "Binding Energy", good: "More negative = stronger binding", desc: "Energy difference between bound and unbound states. Primary ranking metric." },
  { key: "comRMSD", unit: "A", label: "Combined RMSD", good: "Lower is better (<1.0 A ideal)", desc: "How far the placed molecule moved from the original fragment positions. Measures geometric fidelity." },
  { key: "LE", unit: "", label: "Ligand Efficiency", good: "Higher is better (0.3-0.4 drug-like)", desc: "Binding energy normalized by size (-ddG / heavy atoms). Prevents bias toward large molecules." },
  { key: "HA", unit: "", label: "Heavy Atoms", good: "Typically <30 for drug-likeness", desc: "Count of non-hydrogen atoms. Size indicator for drug-likeness assessment." },
  { key: "Const%", unit: "", label: "Constrained Ratio", good: "Higher is better (>50%)", desc: "Fraction of atoms that were constrained to parent fragment positions during placement." },
  { key: "Outcome", unit: "", label: "Outcome Category", good: "'acceptable' is the target", desc: "Classification: acceptable (passed all checks), deviant (moved too much), crashed, timeout, etc." },
];

export function MetricsLegend() {
  const [open, setOpen] = useState(false);

  return (
    <div className="mb-4">
      <button
        onClick={() => setOpen(!open)}
        className="flex items-center gap-1.5 text-[10px] text-slate-400 hover:text-slate-600 transition-colors"
      >
        <i className={`pi ${open ? "pi-chevron-down" : "pi-chevron-right"} text-[8px]`} />
        <i className="pi pi-info-circle text-[10px]" />
        <span className="uppercase tracking-wider font-semibold">Metrics Guide</span>
      </button>
      {open && (
        <div className="mt-2 p-4 rounded-lg bg-slate-50 border border-slate-200">
          <div className="grid gap-3">
            {METRICS.map(m => (
              <div key={m.key} className="flex gap-3">
                <span className="font-mono text-[11px] font-bold text-blue-600 w-14 shrink-0">{m.key}</span>
                <div>
                  <span className="text-xs text-slate-700 font-medium">{m.label}</span>
                  {m.unit && <span className="text-[10px] text-slate-400 ml-1">({m.unit})</span>}
                  <p className="text-[10px] text-slate-400 leading-relaxed">{m.desc}</p>
                  <p className="text-[10px] text-emerald-600">{m.good}</p>
                </div>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
