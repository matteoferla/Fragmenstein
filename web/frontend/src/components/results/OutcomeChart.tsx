"use client";

import { OUTCOME_COLORS, OUTCOME_ORDER } from "@/lib/constants";
import type { ResultRow } from "@/services/types";

interface OutcomeChartProps {
  results: ResultRow[];
}

export function OutcomeChart({ results }: OutcomeChartProps) {
  const counts: Record<string, number> = {};
  for (const r of results) {
    const o = r.outcome || "unknown";
    counts[o] = (counts[o] || 0) + 1;
  }

  const total = results.length;
  if (total === 0) return null;

  return (
    <div className="panel p-5">
      <div className="flex items-center justify-between mb-4">
        <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500">
          Outcome Distribution
        </h3>
        <span className="text-xs font-mono text-slate-400">n={total}</span>
      </div>
      <div className="space-y-2.5">
        {OUTCOME_ORDER.map((outcome) => {
          const count = counts[outcome] || 0;
          if (count === 0) return null;
          const pct = (count / total) * 100;
          const color = OUTCOME_COLORS[outcome] || "#94a3b8";
          return (
            <div key={outcome} className="flex items-center gap-3">
              <span className="text-[10px] w-24 text-right font-mono uppercase tracking-wider text-slate-400">
                {outcome}
              </span>
              <div className="flex-1 h-3 rounded-sm overflow-hidden bg-slate-100">
                <div
                  className="h-full rounded-sm transition-all duration-500"
                  style={{ width: `${pct}%`, background: color }}
                />
              </div>
              <span className="text-[10px] w-16 font-mono" style={{ color }}>
                {count} <span className="text-slate-400">({pct.toFixed(0)}%)</span>
              </span>
            </div>
          );
        })}
      </div>
    </div>
  );
}
