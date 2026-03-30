/** Base URL for API requests. Empty string uses Next.js rewrites proxy. */
export const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || "";

/** Available Victor backend types. */
export const VICTOR_TYPES = ["Wictor", "Victor"] as const;

/** Pipeline step definitions for the step indicator. */
export const STEP_LABELS = [
  { label: "Upload", path: "upload" },
  { label: "Combine", path: "combine" },
  { label: "Similars", path: "similars" },
  { label: "Place", path: "place" },
  { label: "Results", path: "results" },
] as const;

/** Fragmenstein outcome colours — matches LabBench.category_labels order. */
export const OUTCOME_COLORS: Record<string, string> = {
  acceptable: "#22c55e",
  deviant: "#eab308",
  "equally sized": "#f59e0b",
  unstable: "#f97316",
  timeout: "#ef4444",
  "too distant": "#dc2626",
  crashed: "#991b1b",
  unknown: "#94a3b8",
};

/** Display order for outcome categories (best → worst). */
export const OUTCOME_ORDER = [
  "acceptable", "deviant", "equally sized", "unstable",
  "timeout", "too distant", "crashed", "unknown",
] as const;
