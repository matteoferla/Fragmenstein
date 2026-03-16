"use client";

import { useRef, useState, useEffect } from "react";
import { FileUpload, FileUploadHandlerEvent } from "primereact/fileupload";
import { Button } from "primereact/button";
import { InputText } from "primereact/inputtext";
import { Message } from "primereact/message";
import { MolViewer3D } from "@/components/viewer/MolViewer3D";
import { useSessionStore } from "@/stores/sessionStore";
import * as api from "@/services/api";

export function TemplateUpload() {
  const { sessionId, session, loadSession } = useSessionStore();
  const [uploading, setUploading] = useState(false);
  const [message, setMessage] = useState<string | null>(null);
  const [pdbText, setPdbText] = useState<string | null>(null);
  const [removeResidues, setRemoveResidues] = useState("HOH");
  const [cleaning, setCleaning] = useState(false);
  const fileUploadRef = useRef<FileUpload>(null);

  useEffect(() => {
    if (sessionId && session?.template_filename) {
      api.getTemplatePdb(sessionId).then(r => setPdbText(r.pdb)).catch(() => {});
    }
  }, [sessionId, session?.template_filename]);

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

  const handleClean = async () => {
    if (!sessionId || !removeResidues.trim()) return;
    setCleaning(true);
    try {
      const result = await api.cleanTemplate(sessionId, removeResidues);
      setMessage(result.message);
      // Reload PDB for viewer
      const r = await api.getTemplatePdb(sessionId);
      setPdbText(r.pdb);
    } catch (err: unknown) {
      setMessage(err instanceof Error ? err.message : "Clean failed");
    } finally {
      setCleaning(false);
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

      {/* Template preparation: remove unwanted residues */}
      {session?.template_filename && (
        <div className="mt-4 p-3 rounded-lg bg-slate-50 border border-slate-200">
          <div className="text-[10px] font-semibold uppercase tracking-wider text-slate-400 mb-2">
            Template Preparation
          </div>
          <div className="flex items-center gap-2">
            <InputText
              value={removeResidues}
              onChange={(e) => setRemoveResidues(e.target.value.toUpperCase())}
              placeholder="HOH SO4 CL"
              className="flex-1"
            />
            <Button
              label="Remove Residues"
              icon="pi pi-filter"
              size="small"
              severity="secondary"
              onClick={handleClean}
              loading={cleaning}
            />
          </div>
          <div className="text-[9px] text-slate-400 mt-1">
            Remove waters (HOH), ions, and other unwanted residues. Space-separated 3-letter codes.
          </div>
        </div>
      )}

      {/* 3D Protein Viewer */}
      {pdbText && (
        <div className="mt-4">
          <div className="text-[10px] uppercase tracking-wider text-slate-400 mb-2 font-semibold">
            Protein Structure Preview
          </div>
          <MolViewer3D
            models={[{ data: pdbText, format: "pdb" }]}
            height="280px"
          />
        </div>
      )}
    </div>
  );
}
