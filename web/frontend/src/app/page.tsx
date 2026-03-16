"use client";

import { useState, useEffect } from "react";
import { useRouter } from "next/navigation";
import { Button } from "primereact/button";
import { InputText } from "primereact/inputtext";
import { DataTable } from "primereact/datatable";
import { Column } from "primereact/column";
import { useSessionStore } from "@/stores/sessionStore";
import * as api from "@/services/api";
import type { SessionResponse } from "@/services/types";

export default function Home() {
  const router = useRouter();
  const { createSession } = useSessionStore();
  const [sessions, setSessions] = useState<SessionResponse[]>([]);
  const [name, setName] = useState("");
  const [creating, setCreating] = useState(false);

  useEffect(() => {
    api.listSessions().then(setSessions).catch(() => {});
  }, []);

  const handleCreate = async () => {
    setCreating(true);
    try {
      const id = await createSession(name);
      router.push(`/sessions/${id}/upload`);
    } catch {
      // error in store
    } finally {
      setCreating(false);
    }
  };

  const handleResume = (session: SessionResponse) => {
    router.push(`/sessions/${session.id}/upload`);
  };

  return (
    <div className="max-w-4xl mx-auto">
      <div className="mb-10">
        <h1 className="text-3xl font-bold tracking-tight text-slate-900 mb-2">
          Fragmenstein
        </h1>
        <p className="text-sm leading-relaxed max-w-xl text-slate-500">
          Fragment-based drug design pipeline. Upload protein structures and hit
          fragments, combine them into new ligands, search for purchasable
          analogs, and browse ranked results with interactive 3D visualization.
        </p>
      </div>

      {/* New Session */}
      <div className="panel p-6 mb-8">
        <div className="flex items-center gap-2 mb-4">
          <div className="w-6 h-6 rounded-md flex items-center justify-center bg-teal-50 border border-teal-200">
            <i className="pi pi-plus text-xs text-teal-600" />
          </div>
          <h2 className="text-xs font-semibold uppercase tracking-wider text-slate-500">
            New Session
          </h2>
        </div>
        <div className="flex gap-3">
          <InputText
            value={name}
            onChange={(e) => setName(e.target.value)}
            placeholder="Session name (optional)"
            className="flex-1"
          />
          <Button
            label="Initialize"
            icon="pi pi-bolt"
            size="small"
            onClick={handleCreate}
            loading={creating}
          />
        </div>
      </div>

      {/* Recent Sessions */}
      {sessions.length > 0 && (
        <div>
          <div className="flex items-center gap-2 mb-3">
            <div className="w-6 h-6 rounded-md flex items-center justify-center bg-indigo-50 border border-indigo-200">
              <i className="pi pi-history text-xs text-indigo-500" />
            </div>
            <h2 className="text-xs font-semibold uppercase tracking-wider text-slate-500">
              Recent Sessions
            </h2>
          </div>
          <DataTable
            value={sessions}
            onRowClick={(e) => handleResume(e.data as SessionResponse)}
            rowClassName={() => "cursor-pointer"}
            emptyMessage="No sessions yet"
          >
            <Column
              field="name"
              header="Name"
              body={(row: SessionResponse) => (
                <span className="text-slate-800 font-medium">
                  {row.name || (
                    <span className="font-mono text-xs text-slate-400">
                      {row.id.slice(0, 8)}
                    </span>
                  )}
                </span>
              )}
            />
            <Column
              field="status"
              header="Status"
              body={(row: SessionResponse) => (
                <span className="text-xs font-mono px-2 py-0.5 rounded bg-slate-100 text-slate-500 border border-slate-200">
                  {row.status}
                </span>
              )}
            />
            <Column
              field="hit_count"
              header="Hits"
              body={(row: SessionResponse) => (
                <span className="font-mono text-sm text-teal-600 font-semibold">
                  {row.hit_count}
                </span>
              )}
            />
            <Column
              field="created_at"
              header="Created"
              body={(row: SessionResponse) => (
                <span className="text-xs font-mono text-slate-400">
                  {new Date(row.created_at).toLocaleDateString()}
                </span>
              )}
            />
          </DataTable>
        </div>
      )}
    </div>
  );
}
