"use client";

interface MolGrid2DProps {
  smilesList: Array<{ name: string; smiles: string }>;
  columns?: number;
}

export function MolGrid2D({ smilesList, columns = 4 }: MolGrid2DProps) {
  return (
    <div className="grid gap-2" style={{ gridTemplateColumns: `repeat(${columns}, 1fr)` }}>
      {smilesList.map((mol) => (
        <div key={mol.name} className="stat-card text-center">
          <div className="text-xs font-semibold text-slate-700 mb-1">{mol.name}</div>
          <div className="text-[10px] font-mono truncate text-slate-400">{mol.smiles}</div>
        </div>
      ))}
    </div>
  );
}
