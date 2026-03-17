"use client";

import { useEffect, useRef } from "react";

interface SmilesImageProps {
  smiles: string;
  width?: number;
  height?: number;
  className?: string;
  style?: React.CSSProperties;
}

export function SmilesImage({ smiles, width = 200, height = 140, className, style }: SmilesImageProps) {
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!containerRef.current || !smiles) return;
    const el = containerRef.current;
    let cancelled = false;

    // Clear previous children
    while (el.firstChild) el.removeChild(el.firstChild);

    import("smiles-drawer").then((mod) => {
      if (cancelled) return;
      const SmilesDrawer = mod.default;
      const drawer = new SmilesDrawer.SmiDrawer({ width, height });
      const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
      svg.setAttribute("width", String(width));
      svg.setAttribute("height", String(height));
      el.appendChild(svg);
      drawer.draw(smiles, svg, "light", () => {}, () => {
        if (cancelled) return;
        while (el.firstChild) el.removeChild(el.firstChild);
        const span = document.createElement("span");
        span.className = "text-slate-300 text-xs";
        span.textContent = "Invalid SMILES";
        el.appendChild(span);
      });
    });

    return () => { cancelled = true; };
  }, [smiles, width, height]);

  return (
    <div
      ref={containerRef}
      className={className}
      style={{ width, height, display: "flex", alignItems: "center", justifyContent: "center", ...style }}
    />
  );
}
