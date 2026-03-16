"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";

const NAV_ITEMS = [
  { href: "/", label: "Sessions", icon: "pi pi-th-large", exact: true },
];

const TOOL_ITEMS = [
  { path: "monster", label: "Monster", icon: "pi pi-bolt", color: "text-amber-500" },
  { path: "single", label: "Single Victor", icon: "pi pi-wrench", color: "text-violet-500" },
];

export function Sidebar() {
  const pathname = usePathname();

  // Extract session ID from path if in a session
  const sessionMatch = pathname.match(/\/sessions\/([^/]+)/);
  const sessionId = sessionMatch?.[1];

  return (
    <aside className="w-60 flex flex-col bg-white border-r border-slate-200">
      {/* Brand */}
      <div className="p-5 border-b border-slate-100">
        <Link href="/" className="block group">
          <div className="flex items-center gap-2.5">
            <div className="w-8 h-8 rounded-lg flex items-center justify-center text-sm font-bold bg-blue-50 border border-blue-200 text-blue-600">
              F
            </div>
            <div>
              <span className="text-base font-bold tracking-wide text-slate-800 group-hover:text-blue-700 transition-colors">
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
        {NAV_ITEMS.map(item => (
          <Link
            key={item.href}
            href={item.href}
            className={`flex items-center gap-2.5 px-3 py-2 rounded-lg text-sm transition-all ${
              (item.exact ? pathname === item.href : pathname.startsWith(item.href))
                ? "bg-blue-50 text-blue-700 border-l-2 border-blue-500"
                : "text-slate-500 hover:text-slate-700 hover:bg-slate-50"
            }`}
          >
            <i className={`${item.icon} text-xs`} />
            {item.label}
          </Link>
        ))}

        {/* Tools section (only visible in a session) */}
        {sessionId && (
          <>
            <div className="pt-4 pb-1 px-3">
              <span className="text-[9px] font-semibold uppercase tracking-widest text-slate-300">
                Tools
              </span>
            </div>
            {TOOL_ITEMS.map(item => {
              const href = `/sessions/${sessionId}/${item.path}`;
              const isActive = pathname.includes(`/${item.path}`);
              return (
                <Link
                  key={item.path}
                  href={href}
                  className={`flex items-center gap-2.5 px-3 py-2 rounded-lg text-sm transition-all ${
                    isActive
                      ? "bg-slate-50 text-slate-700 border-l-2 border-slate-400"
                      : "text-slate-500 hover:text-slate-700 hover:bg-slate-50"
                  }`}
                >
                  <i className={`${item.icon} text-xs ${item.color}`} />
                  {item.label}
                </Link>
              );
            })}
          </>
        )}
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
