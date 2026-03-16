"use client";

import { useEffect } from "react";
import { useParams } from "next/navigation";
import { StepIndicator } from "@/components/layout/StepIndicator";
import { useSessionStore } from "@/stores/sessionStore";

export default function SessionLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  const params = useParams();
  const sessionId = params.id as string;
  const { loadSession, session } = useSessionStore();

  useEffect(() => {
    if (sessionId) {
      loadSession(sessionId);
    }
  }, [sessionId, loadSession]);

  return (
    <div>
      <div className="flex items-center gap-2 mb-3">
        <div className="w-1.5 h-1.5 rounded-full bg-teal-500" />
        <span className="text-xs font-mono tracking-wider uppercase text-slate-400">
          Session: {session?.name || sessionId.slice(0, 8)}
        </span>
      </div>
      <StepIndicator sessionId={sessionId} />
      {children}
    </div>
  );
}
