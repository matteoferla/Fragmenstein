"use client";

import { DataTable } from "primereact/datatable";
import { Column } from "primereact/column";
import { API_BASE_URL } from "@/lib/constants";

export interface SimilarRow {
  index: number;
  smiles: string | null;
  name: string | null;
  topodist: number | null;
  ecfp4: number | null;
  daylight: number | null;
  dist: number | null;
  mces: number | null;
  query_smiles: string | null;
  alignment: string | null;
  logP: number | null;
  TPSA: number | null;
  similarity_index: number | null;
  [key: string]: unknown;
}

interface SimilarsTableProps {
  results: SimilarRow[];
  onRowSelect?: (row: SimilarRow) => void;
  selectedRow?: SimilarRow | null;
}

function structureImg(row: SimilarRow) {
  if (!row.smiles) return <span className="text-slate-300">-</span>;
  return (
    <img
      src={`${API_BASE_URL}/api/depict?smiles=${encodeURIComponent(row.smiles)}&width=200&height=140`}
      alt={row.smiles}
      className="rounded border border-slate-200"
      style={{ width: 100, height: 70, objectFit: "contain", background: "#fff" }}
      loading="lazy"
    />
  );
}

function numCol(value: number | string | null | undefined, decimals: number = 2) {
  if (value === null || value === undefined) return <span className="text-slate-300">-</span>;
  const num = typeof value === "string" ? parseFloat(value) : value;
  if (isNaN(num)) return <span className="text-slate-300">-</span>;
  return <span className="font-mono text-xs">{num.toFixed(decimals)}</span>;
}

function distBadge(row: SimilarRow) {
  const d = row.topodist;
  if (d === null || d === undefined) return <span className="text-slate-300">-</span>;
  const color = d <= 4 ? "text-emerald-600 bg-emerald-50 border-emerald-200"
    : d <= 8 ? "text-amber-600 bg-amber-50 border-amber-200"
    : "text-red-600 bg-red-50 border-red-200";
  return (
    <span className={`font-mono text-[11px] font-bold px-1.5 py-0.5 rounded border ${color}`}>
      {d}
    </span>
  );
}

export function SimilarsTable({ results, onRowSelect, selectedRow }: SimilarsTableProps) {
  // Detect which columns are present in the data
  const hasCol = (col: string) => results.length > 0 && results.some(r => (r as Record<string, unknown>)[col] != null);

  const hasTopodist = hasCol("topodist");
  const hasTanimotoMerger = hasCol("tanimoto_to_merger");
  const hasEcfp4 = hasCol("ecfp4");
  const hasDaylight = hasCol("daylight");
  const hasQuerySmiles = hasCol("query_smiles");
  const hasMW = hasCol("molecular_weight");
  const hasLogP = hasCol("logP");
  const hasTPSA = hasCol("TPSA");
  const hasSimilarityIndex = hasCol("similarity_index");

  // Pick best sort field
  const defaultSort = hasTanimotoMerger ? "tanimoto_to_merger" : hasSimilarityIndex ? "similarity_index" : hasTopodist ? "topodist" : "name";
  const defaultOrder = hasTanimotoMerger || hasSimilarityIndex ? -1 : hasTopodist ? 1 : 1;

  return (
    <DataTable
      value={results}
      paginator
      rows={15}
      rowsPerPageOptions={[10, 15, 30, 50]}
      selectionMode="single"
      selection={selectedRow}
      onSelectionChange={(e) => onRowSelect?.(e.value as SimilarRow)}
      sortField={defaultSort}
      sortOrder={defaultOrder}
      scrollable
      scrollHeight="700px"
      size="small"
      emptyMessage="No analogs found"
      filterDisplay="row"
    >
      <Column header="Structure" body={structureImg} style={{ minWidth: "120px", padding: "4px 8px" }} />
      <Column field="smiles" header="SMILES" sortable filter filterPlaceholder="Search" style={{ maxWidth: "200px" }}
        body={(r: SimilarRow) => <span className="text-[10px] font-mono truncate block text-slate-500" style={{ maxWidth: "200px" }}>{r.smiles || "-"}</span>}
      />
      <Column field="name" header="Name / ID" sortable filter filterPlaceholder="Search" style={{ minWidth: "110px" }}
        body={(r: SimilarRow) => <span className="text-xs font-mono text-slate-700">{r.name || "-"}</span>}
      />
      {hasTanimotoMerger && (
        <Column field="tanimoto_to_merger" header="Similarity" sortable style={{ minWidth: "90px" }}
          body={(r: SimilarRow) => {
            const v = (r as Record<string, unknown>).tanimoto_to_merger as number | null;
            if (v == null) return <span className="text-slate-300">-</span>;
            const pct = Math.round(v * 100);
            const color = pct >= 70 ? "text-emerald-600" : pct >= 50 ? "text-amber-600" : "text-slate-400";
            return <span className={`font-mono text-xs font-bold ${color}`}>{pct}%</span>;
          }}
        />
      )}
      {hasTopodist && <Column field="topodist" header="Topo Dist" sortable body={distBadge} style={{ minWidth: "90px" }} />}
      {hasEcfp4 && <Column field="ecfp4" header="ECFP4" sortable body={(r: SimilarRow) => numCol(r.ecfp4, 3)} style={{ minWidth: "70px" }} />}
      {hasDaylight && <Column field="daylight" header="Tanimoto" sortable body={(r: SimilarRow) => numCol(r.daylight, 3)} style={{ minWidth: "80px" }} />}
      {hasSimilarityIndex && (
        <Column field="similarity_index" header="Similarity" sortable style={{ minWidth: "90px" }}
          body={(r: SimilarRow) => {
            const v = r.similarity_index;
            if (v == null) return <span className="text-slate-300">-</span>;
            const pct = Math.round(v * 100);
            const color = pct >= 80 ? "text-emerald-600" : pct >= 60 ? "text-amber-600" : "text-slate-400";
            return <span className={`font-mono text-xs font-bold ${color}`}>{pct}%</span>;
          }}
        />
      )}
      {hasMW && (
        <Column field="molecular_weight" header="MW" sortable style={{ minWidth: "70px" }}
          body={(r: SimilarRow) => numCol((r as Record<string, unknown>).molecular_weight as number | null, 1)}
        />
      )}
      {hasLogP && (
        <Column field="logP" header="logP" sortable style={{ minWidth: "70px" }}
          body={(r: SimilarRow) => numCol(r.logP, 2)}
        />
      )}
      {hasTPSA && (
        <Column field="TPSA" header="TPSA" sortable style={{ minWidth: "70px" }}
          body={(r: SimilarRow) => numCol(r.TPSA, 1)}
        />
      )}
      {hasQuerySmiles && (
        <Column field="query_smiles" header="Query (Merger)" style={{ maxWidth: "180px" }}
          body={(r: SimilarRow) => {
            if (!r.query_smiles) return <span className="text-slate-300">-</span>;
            return (
              <div className="flex items-center gap-2">
                <img src={`${API_BASE_URL}/api/depict?smiles=${encodeURIComponent(r.query_smiles)}&width=120&height=80`} alt="query" className="rounded border border-slate-100" style={{ width: 60, height: 40, objectFit: "contain", background: "#fff" }} loading="lazy" />
                <span className="text-[9px] font-mono text-slate-400 truncate" style={{ maxWidth: "100px" }}>{r.query_smiles}</span>
              </div>
            );
          }}
        />
      )}
    </DataTable>
  );
}
