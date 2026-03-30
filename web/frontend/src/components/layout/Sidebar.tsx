"use client";

import { useEffect, useState } from "react";
import Image from "next/image";
import Link from "next/link";
import { usePathname } from "next/navigation";
import * as api from "@/services/api";

const NAV_ITEMS = [
  { href: "/", label: "Sessions", icon: "pi pi-th-large", exact: true },
];

const TOOL_ITEMS = [
  { path: "monster", label: "Monster", icon: "pi pi-bolt", color: "text-amber-500" },
  { path: "single", label: "Single Victor", icon: "pi pi-wrench", color: "text-violet-500" },
];

interface SystemInfo {
  version?: string;
  python?: string;
  platform?: string;
  cores?: number;
  pyrosetta?: boolean;
  pyrosetta_version?: string;
  rdkit?: string;
  gpu?: string;
}

export function Sidebar() {
  const pathname = usePathname();
  const [sysInfo, setSysInfo] = useState<SystemInfo | null>(null);

  useEffect(() => {
    api.getSystemInfo().then((info) => setSysInfo(info as SystemInfo)).catch(() => {});
  }, []);

  const sessionMatch = pathname.match(/\/sessions\/([^/]+)/);
  const sessionId = sessionMatch?.[1];

  return (
    <aside className="w-60 flex flex-col bg-white border-r border-slate-200">
      {/* Brand */}
      <div className="p-5 border-b border-slate-100">
        <Link href="/" className="block group">
          <div className="flex items-center gap-2.5">
            <Image src="/logo.png" alt="Fragmenstein" width={32} height={32} className="rounded-lg" />
            <div>
              <span className="text-base font-bold tracking-wide text-slate-800 group-hover:text-blue-700 transition-colors">
                Fragmenstein
              </span>
              <p className="text-[10px] tracking-widest uppercase text-slate-400">
                Drug Design
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
      <div className="p-4 border-t border-slate-100 space-y-2.5">
        {/* System status */}
        <div className="space-y-1">
          <div className="flex items-center gap-2">
            <div className={`w-1.5 h-1.5 rounded-full ${sysInfo ? "bg-emerald-500" : "bg-slate-300"}`} />
            <span className="text-[10px] tracking-wider uppercase text-slate-400">
              {sysInfo ? "System Ready" : "Connecting..."}
            </span>
            <span className="ml-auto text-[10px] font-mono text-slate-400">v{sysInfo?.version || "..."}</span>
          </div>
          {sysInfo && (
            <div className="grid grid-cols-2 gap-x-2 gap-y-0.5 text-[9px] font-mono text-slate-400 pl-3.5">
              <span>PyRosetta</span>
              <span className={sysInfo.pyrosetta ? "text-emerald-500" : "text-red-400"}>
                {sysInfo.pyrosetta ? "Yes" : "No"}
              </span>
              {sysInfo.rdkit && (
                <>
                  <span>RDKit</span>
                  <span>{sysInfo.rdkit}</span>
                </>
              )}
              <span>CPU Cores</span>
              <span>{sysInfo.cores}</span>
              {sysInfo.gpu && (
                <>
                  <span>GPU</span>
                  <span className="text-emerald-500 truncate">{sysInfo.gpu}</span>
                </>
              )}
              {!sysInfo.gpu && (
                <>
                  <span>GPU</span>
                  <span className="text-slate-300">None</span>
                </>
              )}
            </div>
          )}
        </div>

        {/* Links */}
        <div className="grid grid-cols-2 gap-1">
          <a href="https://github.com/matteoferla/Fragmenstein" target="_blank" rel="noopener noreferrer"
            className="flex items-center justify-center gap-1 text-[9px] text-slate-400 hover:text-slate-600 transition-colors py-1 rounded hover:bg-slate-50">
            <i className="pi pi-github text-[10px]" />
            <span>Core</span>
          </a>
          <a href="https://github.com/sidxz/Fragmenstein" target="_blank" rel="noopener noreferrer"
            className="flex items-center justify-center gap-1 text-[9px] text-slate-400 hover:text-slate-600 transition-colors py-1 rounded hover:bg-slate-50">
            <i className="pi pi-github text-[10px]" />
            <span>Web UI</span>
          </a>
          <a href="https://fragmenstein.readthedocs.io/en/latest/" target="_blank" rel="noopener noreferrer"
            className="flex items-center justify-center gap-1 text-[9px] text-slate-400 hover:text-slate-600 transition-colors py-1 rounded hover:bg-slate-50">
            <i className="pi pi-book text-[10px]" />
            <span>Docs</span>
          </a>
          <a href="https://chemrxiv.org/doi/full/10.26434/chemrxiv-2024-17w01" target="_blank" rel="noopener noreferrer"
            className="flex items-center justify-center gap-1 text-[9px] text-slate-400 hover:text-slate-600 transition-colors py-1 rounded hover:bg-slate-50">
            <i className="pi pi-file text-[10px]" />
            <span>Paper</span>
          </a>
        </div>
      </div>
    </aside>
  );
}
