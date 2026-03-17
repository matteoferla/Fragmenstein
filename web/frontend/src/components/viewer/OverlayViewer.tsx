"use client";

import { useEffect, useRef } from "react";

interface OverlayViewerProps {
  proteinPdb?: string;
  hitMolBlocks?: Array<{ name: string; data: string }>;
  resultMolBlock?: string;
  width?: string;
  height?: string;
}

export function OverlayViewer({
  proteinPdb,
  hitMolBlocks = [],
  resultMolBlock,
  width = "100%",
  height = "500px",
}: OverlayViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!containerRef.current) return;
    if (!proteinPdb && hitMolBlocks.length === 0 && !resultMolBlock) return;

    let viewer: ReturnType<typeof import("3dmol")["createViewer"]> | null = null;

    import("3dmol").then(($3Dmol) => {
      if (!containerRef.current) return;

      while (containerRef.current.firstChild) {
        containerRef.current.removeChild(containerRef.current.firstChild);
      }

      viewer = $3Dmol.createViewer(containerRef.current, {
        backgroundColor: "#f0f4f8",
      });

      if (proteinPdb) {
        const protein = viewer.addModel(proteinPdb, "pdb");
        protein.setStyle({ hetflag: false }, { cartoon: { color: "#94a3b8", opacity: 0.7 } });
        protein.setStyle({ hetflag: true }, { stick: { colorscheme: "whiteCarbon", radius: 0.12, opacity: 0.4 } });
      }

      for (const hit of hitMolBlocks) {
        const m = viewer.addModel(hit.data, "mol");
        m.setStyle({}, { stick: { colorscheme: "greenCarbon", radius: 0.15 } });
      }

      if (resultMolBlock) {
        const result = viewer.addModel(resultMolBlock, "mol");
        result.setStyle({}, { stick: { colorscheme: "cyanCarbon", radius: 0.25 } });
      }

      viewer.zoomTo();
      viewer.render();
    });

    return () => {
      if (viewer) { try { viewer.clear(); } catch {} }
    };
  }, [proteinPdb, hitMolBlocks, resultMolBlock]);

  return (
    <div
      ref={containerRef}
      className="rounded-xl overflow-hidden border border-slate-200 shadow-sm"
      style={{ width, height, position: "relative" }}
    />
  );
}
