"use client";

import { useRef, useState, useEffect, useCallback } from "react";
import { FileUpload, FileUploadHandlerEvent } from "primereact/fileupload";
import { Button } from "primereact/button";
import { InputText } from "primereact/inputtext";
import { InputNumber } from "primereact/inputnumber";
import { Checkbox } from "primereact/checkbox";
import { Message } from "primereact/message";
import { MolViewer3D } from "@/components/viewer/MolViewer3D";
import { JobProgress } from "@/components/jobs/JobProgress";
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

  // PyRosetta detection
  const [hasPyRosetta, setHasPyRosetta] = useState(false);

  // Preparation state
  const [prepJobId, setPrepJobId] = useState<string | null>(null);
  const [preparing, setPreparing] = useState(false);
  const [parameterize, setParameterize] = useState(true);
  const [minimize, setMinimize] = useState(false);
  const [centerResi, setCenterResi] = useState<number | null>(null);
  const [centerChain, setCenterChain] = useState("A");
  const [neighborhoodRadius, setNeighborhoodRadius] = useState(8.0);
  const [cycles, setCycles] = useState(3);

  useEffect(() => {
    api.getSystemInfo().then((info) => {
      setHasPyRosetta(!!info.pyrosetta);
    }).catch(() => {});
  }, []);

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
      const r = await api.getTemplatePdb(sessionId);
      setPdbText(r.pdb);
    } catch (err: unknown) {
      setMessage(err instanceof Error ? err.message : "Clean failed");
    } finally {
      setCleaning(false);
    }
  };

  const handlePrepare = async () => {
    if (!sessionId) return;
    setPreparing(true);
    setMessage(null);
    try {
      const result = await api.prepareTemplate(sessionId, {
        parameterize,
        minimize,
        center_resi: minimize ? centerResi : null,
        center_chain: centerChain,
        neighborhood_radius: neighborhoodRadius,
        cycles,
        remove_residues: removeResidues.trim() ? removeResidues.trim().split(/\s+/) : [],
      });
      setPrepJobId(result.job_id);
    } catch (err: unknown) {
      setMessage(err instanceof Error ? err.message : "Preparation failed");
      setPreparing(false);
    }
  };

  const handlePrepComplete = useCallback(async () => {
    setPreparing(false);
    setMessage("Template prepared successfully");
    if (sessionId) {
      const r = await api.getTemplatePdb(sessionId);
      setPdbText(r.pdb);
    }
  }, [sessionId]);

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

      {/* Template Preparation */}
      {session?.template_filename && (
        <div className="mt-4 p-3 rounded-lg bg-slate-50 border border-slate-200">
          <div className="text-[10px] font-semibold uppercase tracking-wider text-slate-400 mb-3">
            Template Preparation
          </div>

          {/* Remove residues — always available */}
          <div className="flex items-center gap-2 mb-3">
            <InputText
              value={removeResidues}
              onChange={(e) => setRemoveResidues(e.target.value.toUpperCase())}
              placeholder="HOH SO4 CL"
              className="flex-1"
            />
            {!hasPyRosetta && (
              <Button
                label="Remove Residues"
                icon="pi pi-filter"
                size="small"
                severity="secondary"
                onClick={handleClean}
                loading={cleaning}
              />
            )}
          </div>
          <div className="text-[9px] text-slate-400 mb-3">
            Space-separated 3-letter codes for residues to remove (waters, ions, etc.)
          </div>

          {/* PyRosetta options */}
          {hasPyRosetta && (
            <>
              <div className="space-y-3 mb-3">
                {/* Parameterize */}
                <div className="flex items-center gap-2">
                  <Checkbox checked={parameterize} onChange={(e) => setParameterize(!!e.checked)} />
                  <span className="text-xs text-slate-600">Parameterize novel ligands</span>
                  <span className="text-[9px] text-slate-400 ml-auto">Auto-generate topology files for non-standard residues</span>
                </div>

                {/* Minimize */}
                <div className="flex items-center gap-2">
                  <Checkbox checked={minimize} onChange={(e) => setMinimize(!!e.checked)} />
                  <span className="text-xs text-slate-600">Energy minimize around target</span>
                  <span className="text-[9px] text-slate-400 ml-auto">FastRelax around a binding-site residue</span>
                </div>

                {/* Minimization params — show when minimize is checked */}
                {minimize && (
                  <div className="ml-6 grid grid-cols-4 gap-3 p-3 rounded bg-white border border-slate-100">
                    <div className="flex flex-col gap-1">
                      <label className="text-[9px] font-semibold uppercase tracking-wider text-slate-400">Residue #</label>
                      <InputNumber value={centerResi} onValueChange={(e) => setCenterResi(e.value ?? null)} min={1} className="w-full" />
                    </div>
                    <div className="flex flex-col gap-1">
                      <label className="text-[9px] font-semibold uppercase tracking-wider text-slate-400">Chain</label>
                      <InputText value={centerChain} onChange={(e) => setCenterChain(e.target.value.toUpperCase())} maxLength={1} className="w-full" />
                    </div>
                    <div className="flex flex-col gap-1">
                      <label className="text-[9px] font-semibold uppercase tracking-wider text-slate-400">Radius (A)</label>
                      <InputNumber value={neighborhoodRadius} onValueChange={(e) => setNeighborhoodRadius(e.value ?? 8)} min={1} max={30} step={0.5} minFractionDigits={1} className="w-full" />
                    </div>
                    <div className="flex flex-col gap-1">
                      <label className="text-[9px] font-semibold uppercase tracking-wider text-slate-400">Cycles</label>
                      <InputNumber value={cycles} onValueChange={(e) => setCycles(e.value ?? 3)} min={1} max={15} className="w-full" />
                    </div>
                  </div>
                )}
              </div>

              {/* Run button */}
              <Button
                label="Prepare Template"
                icon="pi pi-cog"
                size="small"
                onClick={handlePrepare}
                loading={preparing}
                disabled={preparing || (minimize && !centerResi)}
              />

              {/* Job progress */}
              {prepJobId && preparing && (
                <div className="mt-3">
                  <JobProgress jobId={prepJobId} onComplete={handlePrepComplete} />
                </div>
              )}
            </>
          )}
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
