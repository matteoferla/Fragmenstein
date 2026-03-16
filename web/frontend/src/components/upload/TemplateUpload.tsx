"use client";

import { useRef, useState } from "react";
import { FileUpload, FileUploadHandlerEvent } from "primereact/fileupload";
import { Message } from "primereact/message";
import { useSessionStore } from "@/stores/sessionStore";
import * as api from "@/services/api";

export function TemplateUpload() {
  const { sessionId, session, loadSession } = useSessionStore();
  const [uploading, setUploading] = useState(false);
  const [message, setMessage] = useState<string | null>(null);
  const fileUploadRef = useRef<FileUpload>(null);

  const handleUpload = async (e: FileUploadHandlerEvent) => {
    if (!sessionId || e.files.length === 0) return;
    setUploading(true);
    try {
      const result = await api.uploadTemplate(sessionId, e.files[0]);
      setMessage(result.message);
      await loadSession(sessionId);
      fileUploadRef.current?.clear();
    } catch (err: unknown) {
      setMessage(err instanceof Error ? err.message : "Upload failed");
    } finally {
      setUploading(false);
    }
  };

  return (
    <div className="panel p-5">
      <div className="flex items-center gap-2 mb-3">
        <i className="pi pi-box text-xs text-indigo-500" />
        <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500">
          Template Protein (PDB)
        </h3>
      </div>
      {session?.template_filename && (
        <div className="flex items-center gap-2 mb-3 px-3 py-2 rounded-lg bg-emerald-50 border border-emerald-200">
          <div className="w-1.5 h-1.5 rounded-full bg-emerald-500" />
          <span className="text-xs font-mono text-emerald-700">
            {session.template_filename}
          </span>
        </div>
      )}
      <FileUpload
        ref={fileUploadRef}
        mode="advanced"
        accept=".pdb"
        maxFileSize={100 * 1024 * 1024}
        customUpload
        uploadHandler={handleUpload}
        chooseLabel="Select PDB"
        uploadLabel="Upload"
        auto={false}
        chooseOptions={{ className: "p-button-sm" }}
        uploadOptions={{ className: "p-button-sm" }}
        cancelOptions={{ className: "p-button-sm" }}
        emptyTemplate={
          <p className="text-sm p-4 text-slate-400">
            Drag and drop a PDB file here, or click to browse.
          </p>
        }
      />
      {message && <Message severity="info" text={message} className="mt-3 w-full" />}
    </div>
  );
}
