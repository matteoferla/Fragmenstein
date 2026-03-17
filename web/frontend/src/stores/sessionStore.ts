/** Zustand store for session state. */

import { create } from "zustand";
import type { HitInfo, SessionResponse } from "@/services/types";
import * as api from "@/services/api";

interface SessionState {
  sessionId: string | null;
  session: SessionResponse | null;
  hits: HitInfo[];
  combineJobId: string | null;
  similarsJobId: string | null;
  placeJobId: string | null;
  loading: boolean;
  error: string | null;

  createSession: (name?: string) => Promise<string>;
  loadSession: (id: string) => Promise<void>;
  refreshHits: () => Promise<void>;
  setCombineJobId: (id: string | null) => void;
  setSimilarsJobId: (id: string | null) => void;
  setPlaceJobId: (id: string | null) => void;
  clearError: () => void;
}

export const useSessionStore = create<SessionState>((set, get) => ({
  sessionId: null,
  session: null,
  hits: [],
  combineJobId: null,
  similarsJobId: null,
  placeJobId: null,
  loading: false,
  error: null,

  createSession: async (name = "") => {
    set({ loading: true, error: null });
    try {
      const session = await api.createSession(name);
      set({
        sessionId: session.id,
        session,
        hits: [],
        combineJobId: null,
        similarsJobId: null,
        placeJobId: null,
        loading: false,
      });
      return session.id;
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : "Failed to create session";
      set({ error: msg, loading: false });
      throw e;
    }
  },

  loadSession: async (id: string) => {
    set({ loading: true, error: null });
    try {
      const [session, hitsRes, jobs] = await Promise.all([
        api.getSession(id),
        api.getHits(id),
        api.getSessionJobs(id),
      ]);

      // Restore latest job IDs by type (jobs sorted by created_at DESC)
      const latestJob = (type: string) => {
        const matching = jobs.filter(j => j.type === type);
        return matching.length > 0 ? matching[0].id : null;
      };

      set({
        sessionId: id,
        session,
        hits: hitsRes.hits,
        combineJobId: latestJob("combine"),
        similarsJobId: latestJob("similars"),
        placeJobId: latestJob("place"),
        loading: false,
      });
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : "Failed to load session";
      set({ error: msg, loading: false });
    }
  },

  refreshHits: async () => {
    const { sessionId } = get();
    if (!sessionId) return;
    try {
      const hitsRes = await api.getHits(sessionId);
      set({ hits: hitsRes.hits });
    } catch {
      // silent
    }
  },

  setCombineJobId: (id: string | null) => set({ combineJobId: id }),
  setSimilarsJobId: (id: string | null) => set({ similarsJobId: id }),
  setPlaceJobId: (id: string | null) => set({ placeJobId: id }),
  clearError: () => set({ error: null }),
}));
