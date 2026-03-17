"use client";

import { useEffect, useRef } from "react";

interface ModelSpec {
  data: string;
  format: "pdb" | "mol" | "sdf" | "mol2";
  style?: Record<string, unknown>;
  label?: string;
}

interface MolViewer3DProps {
  models: ModelSpec[];
  width?: string;
  height?: string;
  backgroundColor?: string;
}

export function MolViewer3D({
  models,
  width = "100%",
  height = "400px",
  backgroundColor = "#f0f4f8",
}: MolViewer3DProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<unknown>(null);

  useEffect(() => {
    if (!containerRef.current || models.length === 0) return;

    let viewer: ReturnType<typeof import("3dmol")["createViewer"]> | null = null;

    import("3dmol").then(($3Dmol) => {
      if (!containerRef.current) return;

      while (containerRef.current.firstChild) {
        containerRef.current.removeChild(containerRef.current.firstChild);
      }

      viewer = $3Dmol.createViewer(containerRef.current, { backgroundColor });

      for (const model of models) {
        const m = viewer.addModel(model.data, model.format);
        if (model.style) {
          m.setStyle(model.style);
        } else if (model.format === "pdb") {
          m.setStyle({ hetflag: false }, { cartoon: { color: "#94a3b8", opacity: 0.8 } });
          m.setStyle({ hetflag: true }, { stick: { colorscheme: "whiteCarbon", radius: 0.15 } });
        } else {
          m.setStyle({}, { stick: { radius: 0.2 } });
        }
      }

      viewer.zoomTo();
      viewer.render();
      viewerRef.current = viewer;
    });

    return () => {
      if (viewer) { try { viewer.clear(); } catch {} }
    };
  }, [models, backgroundColor]);

  return (
    <div
      ref={containerRef}
      className="rounded-xl overflow-hidden border border-slate-200 shadow-sm"
      style={{ width, height, position: "relative" }}
    />
  );
}
