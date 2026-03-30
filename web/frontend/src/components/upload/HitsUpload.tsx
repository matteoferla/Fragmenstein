"use client";

import { useRef, useState, useEffect } from "react";
import { FileUpload, FileUploadHandlerEvent } from "primereact/fileupload";
import { DataTable } from "primereact/datatable";
import { Column } from "primereact/column";
import { InputText } from "primereact/inputtext";
import { Checkbox } from "primereact/checkbox";
import { Message } from "primereact/message";
import { MolViewer3D } from "@/components/viewer/MolViewer3D";
import { SmilesImage } from "@/components/viewer/SmilesImage";
import { useSessionStore } from "@/stores/sessionStore";

import * as api from "@/services/api";

export function HitsUpload() {
  const { sessionId, hits, refreshHits } = useSessionStore();
  const [, setUploading] = useState(false);
  const [message, setMessage] = useState<string | null>(null);
  const [hitMolBlocks, setHitMolBlocks] = useState<Array<{ name: string; data: string }>>([]);
  const [proteinPdb, setProteinPdb] = useState<string | null>(null);
  const [ligandResn, setLigandResn] = useState("");
  const [proximityBonding, setProximityBonding] = useState(true);
  const fileUploadRef = useRef<FileUpload>(null);

  // Load hit mol blocks + protein for 3D overlay
  useEffect(() => {
    if (!sessionId || hits.length === 0) { setHitMolBlocks([]); return; }
    api.getAllHitMolBlocks(sessionId)
      .then(r => setHitMolBlocks(r.hits.map(h => ({ name: h.name, data: h.mol_block }))))
      .catch(() => {});
    api.getTemplatePdb(sessionId)
      .then(r => setProteinPdb(r.pdb))
      .catch(() => setProteinPdb(null));
  }, [sessionId, hits.length]);

  const handleUpload = async (e: FileUploadHandlerEvent) => {
    if (!sessionId || e.files.length === 0) return;
    setUploading(true);
    try {
      const result = await api.uploadHits(sessionId, e.files, ligandResn || undefined, proximityBonding);
      setMessage(result.message);
      await refreshHits();
      fileUploadRef.current?.clear();
    } catch (err: unknown) {
      setMessage(err instanceof Error ? err.message : "Upload failed");
    } finally {
      setUploading(false);
    }
  };

  // Build 3D viewer models: protein (if available) + all hits
  const viewerModels = [
    ...(proteinPdb ? [{ data: proteinPdb, format: "pdb" as const }] : []),
    ...hitMolBlocks.map(h => ({ data: h.data, format: "mol" as const })),
  ];

  return (
    <div className="panel p-5">
      <div className="flex items-center gap-2 mb-3">
        <i className="pi pi-sitemap text-xs text-blue-600" />
        <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500">
          Hit Compounds (SDF/MOL/PDB)
        </h3>
      </div>
      {/* Ligand residue name — needed when hits are in PDB format */}
      <div className="flex items-center gap-3 mb-3 p-3 rounded-lg bg-slate-50 border border-slate-200">
        <div className="flex-1">
          <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400 block mb-1">
            PDB Ligand Residue Name
          </label>
          <InputText
            value={ligandResn}
            onChange={(e) => setLigandResn(e.target.value.toUpperCase())}
            placeholder="e.g. LIG, UNL, DRG"
            className="w-full"
            maxLength={3}
          />
        </div>
        <div className="flex flex-col gap-2">
          <label className="flex items-center gap-2 cursor-pointer">
            <Checkbox checked={proximityBonding} onChange={(e) => setProximityBonding(e.checked ?? true)} />
            <span className="text-[10px] text-slate-500">Proximity Bonding</span>
          </label>
          <div className="text-[10px] text-slate-400 max-w-[200px]">
            For PDB files: residue name extracts the ligand. Enable proximity bonding if the PDB lacks CONECT records.
          </div>
        </div>
      </div>
      <div className="text-[10px] text-slate-400 mb-3 px-1">
        <i className="pi pi-info-circle mr-1" />
        For PDB hits, you can also include a <strong>.smi</strong> file (tab-separated SMILES and name) to correct bond orders during ligand extraction.
      </div>

      <FileUpload
        ref={fileUploadRef}
        mode="advanced"
        accept=".sdf,.mol,.mol2,.pdb,.smi"
        multiple
        maxFileSize={100 * 1024 * 1024}
        customUpload
        uploadHandler={handleUpload}
        chooseLabel="Select Hits"
        uploadLabel="Upload"
        auto={false}
        chooseOptions={{ className: "p-button-sm" }}
        uploadOptions={{ className: "p-button-sm" }}
        cancelOptions={{ className: "p-button-sm" }}
        emptyTemplate={
          <p className="text-sm p-4 text-slate-400">
            Drag and drop hit molecule files here. Supports SDF, MOL, MOL2, PDB, SMI.
          </p>
        }
      />
      {message && <Message severity="info" text={message} className="mt-3 w-full" />}

      {hits.length > 0 && (
        <div className="mt-5">
          <div className="flex items-center gap-2 mb-2">
            <span className="font-mono text-sm font-bold text-blue-600">{hits.length}</span>
            <span className="text-xs uppercase tracking-wider text-slate-400">hit(s) loaded</span>
          </div>

          {/* 2D Structures grid */}
          <div className="grid grid-cols-4 gap-2 mb-4">
            {hits.map((h, i) => (
              <div key={`${h.name}-${i}`} className="bg-white rounded-lg border border-slate-200 p-2 text-center">
                {h.smiles ? (
                  <SmilesImage smiles={h.smiles} width={200} height={140} className="mx-auto rounded" style={{ width: "100%", height: 70, objectFit: "contain" }} />
                ) : (
                  <div className="h-[70px] flex items-center justify-center text-slate-300 text-xs">No SMILES</div>
                )}
                <div className="text-[10px] font-mono text-slate-600 mt-1 truncate">{h.name}</div>
                <div className="text-[9px] text-slate-400">{h.num_atoms} HA</div>
              </div>
            ))}
          </div>

          {/* 3D overlay: protein + all hits */}
          {viewerModels.length > 0 && (
            <div>
              <div className="text-[10px] uppercase tracking-wider text-slate-400 mb-2 font-semibold">
                3D Overlay — {proteinPdb ? "Protein + " : ""}{hits.length} Hit(s)
              </div>
              <MolViewer3D models={viewerModels} height="300px" />
            </div>
          )}

          {/* Table */}
          <div className="mt-4">
            <DataTable value={hits} size="small" scrollable scrollHeight="200px">
              <Column field="name" header="Name" sortable />
              <Column field="smiles" header="SMILES" style={{ maxWidth: "300px" }}
                body={(row) => (
                  <span className="text-xs font-mono truncate block text-slate-400" style={{ maxWidth: "300px" }}>
                    {row.smiles || "-"}
                  </span>
                )}
              />
              <Column field="num_atoms" header="Heavy Atoms" sortable
                body={(row) => <span className="font-mono text-xs text-blue-600">{row.num_atoms}</span>}
              />
              <Column field="filename" header="Source File" />
            </DataTable>
          </div>
        </div>
      )}
    </div>
  );
}
