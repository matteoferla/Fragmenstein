"use client";

import { useState, useEffect } from "react";
import { useParams } from "next/navigation";
import { TabView, TabPanel } from "primereact/tabview";
import { ResultsTable } from "@/components/results/ResultsTable";
import { OutcomeChart } from "@/components/results/OutcomeChart";
import { ResultDetail } from "@/components/results/ResultDetail";
import { DownloadPanel } from "@/components/results/DownloadPanel";
import { useSessionStore } from "@/stores/sessionStore";
import * as api from "@/services/api";
import type { ResultRow } from "@/services/types";

export default function ResultsPage() {
  const params = useParams();
  const sessionId = params.id as string;
  const { combineJobId, placeJobId } = useSessionStore();

  const [combineResults, setCombineResults] = useState<ResultRow[]>([]);
  const [placeResults, setPlaceResults] = useState<ResultRow[]>([]);
  const [selectedRow, setSelectedRow] = useState<ResultRow | null>(null);
  const [activeJobId, setActiveJobId] = useState<string | null>(null);
  const [activeTab, setActiveTab] = useState(0);

  useEffect(() => {
    if (combineJobId) {
      api.getJobStatus(combineJobId).then((job) => {
        if (job.status === "completed") {
          api.getJobResults(combineJobId).then((res) => setCombineResults(res.results)).catch(() => {});
        }
      }).catch(() => {});
    }
    if (placeJobId) {
      api.getJobStatus(placeJobId).then((job) => {
        if (job.status === "completed") {
          api.getJobResults(placeJobId).then((res) => setPlaceResults(res.results)).catch(() => {});
        }
      }).catch(() => {});
    }
  }, [combineJobId, placeJobId]);

  useEffect(() => { setActiveJobId(combineJobId); }, [combineJobId]);

  const handleTabChange = (index: number) => {
    setActiveTab(index);
    setSelectedRow(null);
    setActiveJobId(index === 0 ? combineJobId : placeJobId);
  };

  const renderPanel = (results: ResultRow[], jobId: string | null, emptyMsg: string) => {
    if (results.length === 0) {
      return (
        <div className="panel p-8 text-center mt-4">
          <i className="pi pi-inbox text-xl mb-2 text-slate-300" />
          <p className="text-sm text-slate-400">{emptyMsg}</p>
        </div>
      );
    }
    return (
      <>
        <div className="mb-4 flex items-center justify-between mt-4">
          <div className="flex items-center gap-2">
            <span className="font-mono text-lg font-bold text-blue-600">{results.length}</span>
            <span className="text-xs uppercase tracking-wider text-slate-400">results</span>
          </div>
          {jobId && <DownloadPanel jobId={jobId} />}
        </div>
        <OutcomeChart results={results} />
        <div className="grid grid-cols-3 gap-5 mt-5">
          <div className="col-span-2">
            <ResultsTable results={results} onRowSelect={(row) => { setSelectedRow(row); setActiveJobId(jobId); }} selectedRow={selectedRow} />
          </div>
          <div>
            {selectedRow && activeJobId === jobId ? (
              <ResultDetail result={selectedRow} jobId={jobId!} sessionId={sessionId} />
            ) : (
              <div className="panel p-8 flex flex-col items-center justify-center text-center" style={{ minHeight: "300px" }}>
                <i className="pi pi-eye text-2xl mb-3 text-slate-300" />
                <p className="text-xs text-slate-400">Select a result to inspect</p>
              </div>
            )}
          </div>
        </div>
      </>
    );
  };

  return (
    <div className="max-w-6xl">
      <div className="flex items-center gap-3 mb-6">
        <div className="w-8 h-8 rounded-lg flex items-center justify-center bg-blue-50 border border-blue-200 text-blue-600">
          <i className="pi pi-chart-bar text-sm" />
        </div>
        <div>
          <h2 className="text-lg font-bold text-slate-800">Results Browser</h2>
          <p className="text-xs text-slate-400">Explore, compare, and export all pipeline results</p>
        </div>
      </div>

      <TabView activeIndex={activeTab} onTabChange={(e) => handleTabChange(e.index)}>
        <TabPanel header={`Combine (${combineResults.length})`}>
          {renderPanel(combineResults, combineJobId, "No combine results available.")}
        </TabPanel>
        <TabPanel header={`Placement (${placeResults.length})`}>
          {renderPanel(placeResults, placeJobId, "No placement results available.")}
        </TabPanel>
      </TabView>
    </div>
  );
}
