"use client";

import { Button } from "primereact/button";
import { getResultsDownloadUrl } from "@/services/api";

interface DownloadPanelProps {
  jobId: string;
}

export function DownloadPanel({ jobId }: DownloadPanelProps) {
  return (
    <div className="flex gap-2">
      <a href={getResultsDownloadUrl(jobId, "csv")} download>
        <Button
          label="Export CSV"
          icon="pi pi-download"
          severity="secondary"
          size="small"
        />
      </a>
    </div>
  );
}
