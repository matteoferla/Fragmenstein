"use client";

import { Button } from "primereact/button";
import { getResultsDownloadUrl } from "@/services/api";

interface DownloadPanelProps {
  jobId: string;
  showSdf?: boolean;
  showPdb?: boolean;
}

export function DownloadPanel({ jobId, showSdf = true, showPdb = true }: DownloadPanelProps) {
  return (
    <div className="flex gap-2">
      <a href={getResultsDownloadUrl(jobId, "csv")} download>
        <Button label="CSV" icon="pi pi-download" severity="secondary" size="small" />
      </a>
      {showSdf && (
        <a href={getResultsDownloadUrl(jobId, "sdf")} download>
          <Button label="SDF" icon="pi pi-download" severity="secondary" size="small" />
        </a>
      )}
      {showPdb && (
        <a href={getResultsDownloadUrl(jobId, "pdb")} download>
          <Button label="PDB (ZIP)" icon="pi pi-download" severity="secondary" size="small" />
        </a>
      )}
    </div>
  );
}
