"use client";

interface JobStatusBadgeProps {
  status: string;
}

const STATUS_STYLES: Record<string, string> = {
  pending: "bg-slate-100 text-slate-500 border-slate-200",
  running: "bg-teal-50 text-teal-700 border-teal-200",
  completed: "bg-emerald-50 text-emerald-700 border-emerald-200",
  failed: "bg-red-50 text-red-700 border-red-200",
  cancelled: "bg-amber-50 text-amber-700 border-amber-200",
};

export function JobStatusBadge({ status }: JobStatusBadgeProps) {
  const style = STATUS_STYLES[status] || STATUS_STYLES.pending;
  return (
    <span className={`px-2.5 py-1 rounded-md text-[10px] font-bold uppercase tracking-wider font-mono border ${style}`}>
      {status}
    </span>
  );
}
