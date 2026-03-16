/** SSE EventSource hook for job progress. */

import { useEffect, useRef, useState, useCallback } from "react";
import { API_BASE_URL } from "@/lib/constants";
import type { JobProgressEvent } from "@/services/types";

interface UseJobProgressReturn {
  progress: number;
  status: string;
  message: string;
  isComplete: boolean;
  isFailed: boolean;
}

export function useJobProgress(jobId: string | null): UseJobProgressReturn {
  const [progress, setProgress] = useState(0);
  const [status, setStatus] = useState("pending");
  const [message, setMessage] = useState("");
  const sourceRef = useRef<EventSource | null>(null);

  useEffect(() => {
    if (!jobId) return;

    const source = new EventSource(`${API_BASE_URL}/api/jobs/${jobId}/stream`);
    sourceRef.current = source;

    source.onmessage = (event) => {
      try {
        const data: JobProgressEvent = JSON.parse(event.data);
        setProgress(data.progress);
        setStatus(data.status);
        setMessage(data.message);

        if (data.status === "completed" || data.status === "failed") {
          source.close();
        }
      } catch {
        // ignore parse errors
      }
    };

    source.onerror = () => {
      source.close();
    };

    return () => {
      source.close();
      sourceRef.current = null;
    };
  }, [jobId]);

  return {
    progress,
    status,
    message,
    isComplete: status === "completed",
    isFailed: status === "failed",
  };
}
