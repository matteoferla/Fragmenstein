"use client";

import { useState } from "react";
import { useParams, useRouter } from "next/navigation";
import { Button } from "primereact/button";
import { TemplateUpload } from "@/components/upload/TemplateUpload";
import { HitsUpload } from "@/components/upload/HitsUpload";
import { useSessionStore } from "@/stores/sessionStore";
import * as api from "@/services/api";

export default function UploadPage() {
  const params = useParams();
  const router = useRouter();
  const { session, hits, loadSession } = useSessionStore();
  const sessionId = params.id as string;
  const [loadingDemo, setLoadingDemo] = useState<string | null>(null);

  const canProceed = session?.template_filename && hits.length >= 2;

  const handleLoadDemo = async (dataset: string) => {
    setLoadingDemo(dataset);
    try {
      await api.loadDemo(sessionId, dataset, 5);
      await loadSession(sessionId);
    } catch {
      // error
    } finally {
      setLoadingDemo(null);
    }
  };

  return (
    <div className="max-w-4xl">
      <div className="flex items-center gap-3 mb-6">
        <div className="w-8 h-8 rounded-lg flex items-center justify-center bg-blue-50 border border-blue-200 text-blue-600">
          <i className="pi pi-upload text-sm" />
        </div>
        <div>
          <h2 className="text-lg font-bold text-slate-800">Upload Structure Data</h2>
          <p className="text-xs text-slate-400">Provide a template protein and hit fragment molecules</p>
        </div>
      </div>

      {/* Demo Quick Start */}
      <div className="panel p-5 mb-6">
        <div className="flex items-center gap-2 mb-3">
          <i className="pi pi-bolt text-xs text-amber-500" />
          <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500">
            Quick Start — Load Demo Dataset
          </h3>
        </div>
        <p className="text-xs text-slate-400 mb-3">
          Load a pre-built template + 5 hit fragments from a published dataset to try the pipeline immediately.
        </p>
        <div className="flex gap-2">
          <Button
            label="MPro (COVID)"
            icon="pi pi-database"
            size="small"
            severity="secondary"
            onClick={() => handleLoadDemo("mpro")}
            loading={loadingDemo === "mpro"}
            disabled={!!loadingDemo}
          />
          <Button
            label="Mac1 (SAR-CoV-2)"
            icon="pi pi-database"
            size="small"
            severity="secondary"
            onClick={() => handleLoadDemo("mac1")}
            loading={loadingDemo === "mac1"}
            disabled={!!loadingDemo}
          />
        </div>
      </div>

      <div className="space-y-6">
        <TemplateUpload />
        <HitsUpload />
      </div>

      <div className="mt-8 flex justify-end">
        <Button
          label="Proceed to Combine"
          icon="pi pi-arrow-right"
          iconPos="right"
          size="small"
          disabled={!canProceed}
          onClick={() => router.push(`/sessions/${sessionId}/combine`)}
          tooltip={!canProceed ? "Upload a template PDB and at least 2 hit molecules" : undefined}
        />
      </div>
    </div>
  );
}
