"use client";

import { Button } from "primereact/button";
import { ProgressBar } from "primereact/progressbar";
import { useJobProgress } from "@/hooks/useJobProgress";
import { JobStatusBadge } from "./JobStatusBadge";
import * as api from "@/services/api";

interface JobProgressProps {
  jobId: string | null;
  onComplete?: () => void;
  onCancel?: () => void;
  onRerun?: () => void;
}

export function JobProgress({ jobId, onComplete, onCancel, onRerun }: JobProgressProps) {
  const { status, message, isComplete, isFailed } =
    useJobProgress(jobId);

  if (isComplete && onComplete) {
    onComplete();
  }

  if (!jobId) return null;

  const isCancelled = status === "cancelled";
  const isRunning = status === "running" || status === "pending";
  const isDone = isComplete || isFailed || isCancelled;

  const handleCancel = async () => {
    try {
      await api.cancelJob(jobId);
      onCancel?.();
    } catch {
      // ignore
    }
  };

  return (
    <div className={`panel p-5 ${isRunning ? "pulse-border" : ""}`}>
      <div className="flex items-center justify-between mb-3">
        <div className="flex items-center gap-2">
          {isRunning && (
            <div className="w-2 h-2 rounded-full bg-blue-500 animate-pulse" />
          )}
          <span className="text-xs font-semibold uppercase tracking-wider text-slate-500">
            Job Progress
          </span>
        </div>
        <div className="flex items-center gap-2">
          <JobStatusBadge status={status} />
          {isRunning && (
            <Button
              icon="pi pi-times"
              label="Cancel"
              size="small"
              severity="danger"
              text
              className="!text-[10px] !py-1 !px-2"
              onClick={handleCancel}
            />
          )}
          {isDone && onRerun && (
            <Button
              icon="pi pi-refresh"
              label="Rerun"
              size="small"
              severity="secondary"
              text
              className="!text-[10px] !py-1 !px-2"
              onClick={onRerun}
            />
          )}
        </div>
      </div>

      <div className="mb-2">
        <div className="flex items-center justify-between mb-1">
          <span className="text-xs font-mono text-slate-400">
            Processing...
          </span>
        </div>
        <ProgressBar mode="indeterminate" style={{ height: "6px" }} />
      </div>

      {message && (
        <p className="text-xs font-mono mt-2 text-slate-400">{message}</p>
      )}
      {isFailed && (
        <p className="text-xs mt-2 text-red-600">
          Job failed. {onRerun ? "Click Rerun to try again." : "Check error details."}
        </p>
      )}
      {isCancelled && (
        <p className="text-xs mt-2 text-amber-600">
          Job was cancelled. {onRerun ? "Click Rerun to start a new run." : ""}
        </p>
      )}
    </div>
  );
}
