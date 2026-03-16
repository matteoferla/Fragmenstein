"use client";

import { useParams, useRouter } from "next/navigation";
import { Button } from "primereact/button";
import { TemplateUpload } from "@/components/upload/TemplateUpload";
import { HitsUpload } from "@/components/upload/HitsUpload";
import { useSessionStore } from "@/stores/sessionStore";

export default function UploadPage() {
  const params = useParams();
  const router = useRouter();
  const { session, hits } = useSessionStore();
  const sessionId = params.id as string;

  const canProceed = session?.template_filename && hits.length >= 2;

  return (
    <div className="max-w-4xl">
      <div className="flex items-center gap-3 mb-6">
        <div className="w-8 h-8 rounded-lg flex items-center justify-center bg-teal-50 border border-teal-200 text-teal-600">
          <i className="pi pi-upload text-sm" />
        </div>
        <div>
          <h2 className="text-lg font-bold text-slate-800">Upload Structure Data</h2>
          <p className="text-xs text-slate-400">Provide a template protein and hit fragment molecules</p>
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
