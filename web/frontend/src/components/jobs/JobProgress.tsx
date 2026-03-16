"use client";

import { ProgressBar } from "primereact/progressbar";
import { useJobProgress } from "@/hooks/useJobProgress";
import { JobStatusBadge } from "./JobStatusBadge";

interface JobProgressProps {
  jobId: string | null;
  onComplete?: () => void;
}

export function JobProgress({ jobId, onComplete }: JobProgressProps) {
  const { progress, status, message, isComplete, isFailed } =
    useJobProgress(jobId);

  if (isComplete && onComplete) {
    onComplete();
  }

  if (!jobId) return null;

  return (
    <div className={`panel p-5 ${status === "running" ? "pulse-border" : ""}`}>
      <div className="flex items-center justify-between mb-3">
        <div className="flex items-center gap-2">
          {status === "running" && (
            <div className="w-2 h-2 rounded-full bg-teal-500 animate-pulse" />
          )}
          <span className="text-xs font-semibold uppercase tracking-wider text-slate-500">
            Job Progress
          </span>
        </div>
        <JobStatusBadge status={status} />
      </div>

      <div className="mb-2">
        <div className="flex items-center justify-between mb-1">
          <span className="text-xs font-mono text-slate-400">
            {Math.round(progress * 100)}%
          </span>
        </div>
        <ProgressBar value={Math.round(progress * 100)} showValue={false} />
      </div>

      {message && (
        <p className="text-xs font-mono mt-2 text-slate-400">{message}</p>
      )}
      {isFailed && (
        <p className="text-xs mt-2 text-red-600">
          Job failed. Check the error details.
        </p>
      )}
    </div>
  );
}
