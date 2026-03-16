"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";

export function Sidebar() {
  const pathname = usePathname();

  return (
    <aside className="w-60 flex flex-col bg-white border-r border-slate-200">
      {/* Brand */}
      <div className="p-5 border-b border-slate-100">
        <Link href="/" className="block group">
          <div className="flex items-center gap-2.5">
            <div className="w-8 h-8 rounded-lg flex items-center justify-center text-sm font-bold bg-teal-50 border border-teal-200 text-teal-600">
              F
            </div>
            <div>
              <span className="text-base font-bold tracking-wide text-slate-800 group-hover:text-teal-700 transition-colors">
                Fragmenstein
              </span>
              <p className="text-[10px] tracking-widest uppercase text-slate-400">
                Drug Design Lab
              </p>
            </div>
          </div>
        </Link>
      </div>

      {/* Nav */}
      <nav className="flex-1 p-3 space-y-1">
        <Link
          href="/"
          className={`flex items-center gap-2.5 px-3 py-2 rounded-lg text-sm transition-all ${
            pathname === "/"
              ? "bg-teal-50 text-teal-700 border-l-2 border-teal-500"
              : "text-slate-500 hover:text-slate-700 hover:bg-slate-50"
          }`}
        >
          <i className="pi pi-th-large text-xs" />
          Sessions
        </Link>
      </nav>

      {/* Footer */}
      <div className="p-4 flex items-center gap-2 border-t border-slate-100 text-slate-400">
        <div className="w-1.5 h-1.5 rounded-full bg-emerald-500" />
        <span className="text-[10px] tracking-wider uppercase">System Ready</span>
        <span className="ml-auto text-[10px] font-mono">v0.1.0</span>
      </div>
    </aside>
  );
}
