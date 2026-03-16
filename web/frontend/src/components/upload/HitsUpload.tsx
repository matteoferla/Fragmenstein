"use client";

import { useRef, useState } from "react";
import { FileUpload, FileUploadHandlerEvent } from "primereact/fileupload";
import { DataTable } from "primereact/datatable";
import { Column } from "primereact/column";
import { Message } from "primereact/message";
import { useSessionStore } from "@/stores/sessionStore";
import * as api from "@/services/api";

export function HitsUpload() {
  const { sessionId, hits, refreshHits } = useSessionStore();
  const [uploading, setUploading] = useState(false);
  const [message, setMessage] = useState<string | null>(null);
  const fileUploadRef = useRef<FileUpload>(null);

  const handleUpload = async (e: FileUploadHandlerEvent) => {
    if (!sessionId || e.files.length === 0) return;
    setUploading(true);
    try {
      const result = await api.uploadHits(sessionId, e.files);
      setMessage(result.message);
      await refreshHits();
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
        <i className="pi pi-sitemap text-xs text-teal-600" />
        <h3 className="text-xs font-semibold uppercase tracking-wider text-slate-500">
          Hit Compounds (SDF/MOL/PDB)
        </h3>
      </div>
      <FileUpload
        ref={fileUploadRef}
        mode="advanced"
        accept=".sdf,.mol,.mol2,.pdb"
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
            Drag and drop hit molecule files here. Supports SDF, MOL, MOL2, PDB.
          </p>
        }
      />
      {message && <Message severity="info" text={message} className="mt-3 w-full" />}

      {hits.length > 0 && (
        <div className="mt-5">
          <div className="flex items-center gap-2 mb-2">
            <span className="font-mono text-sm font-bold text-teal-600">{hits.length}</span>
            <span className="text-xs uppercase tracking-wider text-slate-400">hit(s) loaded</span>
          </div>
          <DataTable value={hits} size="small" scrollable scrollHeight="300px">
            <Column field="name" header="Name" sortable />
            <Column field="smiles" header="SMILES" style={{ maxWidth: "300px" }}
              body={(row) => (
                <span className="text-xs font-mono truncate block text-slate-400" style={{ maxWidth: "300px" }}>
                  {row.smiles || "-"}
                </span>
              )}
            />
            <Column field="num_atoms" header="Heavy Atoms" sortable
              body={(row) => <span className="font-mono text-xs text-teal-600">{row.num_atoms}</span>}
            />
            <Column field="filename" header="Source File" />
          </DataTable>
        </div>
      )}
    </div>
  );
}
