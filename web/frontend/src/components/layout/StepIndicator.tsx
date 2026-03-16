"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";
import { STEP_LABELS } from "@/lib/constants";

interface StepIndicatorProps {
  sessionId: string;
}

export function StepIndicator({ sessionId }: StepIndicatorProps) {
  const pathname = usePathname();

  const currentStep = STEP_LABELS.findIndex((s) =>
    pathname.includes(`/${s.path}`)
  );

  return (
    <div className="flex items-center gap-0 mb-8 p-1 rounded-xl bg-white border border-slate-200">
      {STEP_LABELS.map((step, i) => {
        const isActive = pathname.includes(`/${step.path}`);
        const isPast = i < currentStep;
        const href = `/sessions/${sessionId}/${step.path}`;

        return (
          <div key={step.path} className="flex items-center flex-1">
            {i > 0 && (
              <div className={`w-px h-5 ${isPast ? "bg-teal-300" : "bg-slate-200"}`} />
            )}
            <Link
              href={href}
              className={`flex items-center justify-center gap-2 px-4 py-2 rounded-lg text-xs font-medium transition-all flex-1 text-center ${
                isActive
                  ? "bg-teal-50 text-teal-700 border-b-2 border-teal-500"
                  : isPast
                  ? "text-slate-700"
                  : "text-slate-400 hover:text-slate-600"
              }`}
            >
              <span
                className={`w-5 h-5 rounded-full flex items-center justify-center text-[10px] font-bold ${
                  isActive
                    ? "bg-teal-600 text-white"
                    : isPast
                    ? "bg-teal-100 text-teal-700"
                    : "bg-slate-100 text-slate-400"
                }`}
              >
                {isPast ? "\u2713" : i + 1}
              </span>
              {step.label}
            </Link>
          </div>
        );
      })}
    </div>
  );
}
