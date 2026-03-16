"use client";

import { DataTable } from "primereact/datatable";
import { Column } from "primereact/column";
import { OUTCOME_COLORS } from "@/lib/constants";
import type { ResultRow } from "@/services/types";

interface ResultsTableProps {
  results: ResultRow[];
  onRowSelect?: (row: ResultRow) => void;
  selectedRow?: ResultRow | null;
}

function outcomeBadge(row: ResultRow) {
  const color = OUTCOME_COLORS[row.outcome || ""] || "#94a3b8";
  return (
    <span
      className="text-[10px] font-bold uppercase tracking-wider font-mono px-2 py-0.5 rounded"
      style={{ background: `${color}15`, color, border: `1px solid ${color}30` }}
    >
      {row.outcome || "unknown"}
    </span>
  );
}

function numCol(value: number | null, decimals: number = 2) {
  if (value === null || value === undefined) return (
    <span className="text-slate-300">-</span>
  );
  return <span className="font-mono text-xs">{value.toFixed(decimals)}</span>;
}

export function ResultsTable({ results, onRowSelect, selectedRow }: ResultsTableProps) {
  return (
    <DataTable
      value={results}
      paginator
      rows={20}
      rowsPerPageOptions={[10, 20, 50, 100]}
      selectionMode="single"
      selection={selectedRow}
      onSelectionChange={(e) => onRowSelect?.(e.value as ResultRow)}
      sortField="ddG"
      sortOrder={1}
      scrollable
      scrollHeight="600px"
      size="small"
      emptyMessage="No results"
      filterDisplay="row"
    >
      <Column field="name" header="Name" sortable filter filterPlaceholder="Search" style={{ minWidth: "120px" }}
        body={(r: ResultRow) => <span className="font-mono text-xs text-slate-700">{r.name}</span>}
      />
      <Column field="outcome" header="Outcome" body={outcomeBadge} sortable filter filterPlaceholder="Filter" style={{ minWidth: "110px" }} />
      <Column field="ddG" header="ddG" sortable body={(r: ResultRow) => numCol(r.ddG)} style={{ minWidth: "80px" }} />
      <Column field="comRMSD" header="comRMSD" sortable body={(r: ResultRow) => numCol(r.comRMSD)} style={{ minWidth: "90px" }} />
      <Column field="LE" header="LE" sortable body={(r: ResultRow) => numCol(r.LE, 3)} style={{ minWidth: "70px" }} />
      <Column field="N_HA" header="HA" sortable style={{ minWidth: "60px" }}
        body={(r: ResultRow) => <span className="font-mono text-xs text-blue-600">{r.N_HA ?? "-"}</span>}
      />
      <Column field="smiles" header="SMILES" style={{ maxWidth: "200px" }}
        body={(r: ResultRow) => (
          <span className="text-[11px] font-mono truncate block text-slate-400" style={{ maxWidth: "200px" }}>
            {r.smiles || "-"}
          </span>
        )}
      />
      <Column field="runtime" header="Time (s)" sortable body={(r: ResultRow) => numCol(r.runtime, 1)} style={{ minWidth: "80px" }} />
    </DataTable>
  );
}
