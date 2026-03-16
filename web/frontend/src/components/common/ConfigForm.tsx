"use client";

import { InputText } from "primereact/inputtext";
import { InputNumber } from "primereact/inputnumber";
import { Dropdown } from "primereact/dropdown";
import { Checkbox } from "primereact/checkbox";
import { VICTOR_TYPES } from "@/lib/constants";

interface ConfigField {
  key: string;
  label: string;
  type: "text" | "number" | "select" | "checkbox";
  options?: string[];
  min?: number;
  max?: number;
  step?: number;
}

interface ConfigFormProps {
  fields: ConfigField[];
  values: Record<string, unknown>;
  onChange: (key: string, value: unknown) => void;
}

export function ConfigForm({ fields, values, onChange }: ConfigFormProps) {
  return (
    <div className="grid grid-cols-2 gap-4">
      {fields.map((field) => (
        <div key={field.key} className="flex flex-col gap-1.5">
          <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">
            {field.label}
          </label>
          {field.type === "text" && (
            <InputText
              value={(values[field.key] as string) || ""}
              onChange={(e) => onChange(field.key, e.target.value)}
              className="w-full"
            />
          )}
          {field.type === "number" && (
            <InputNumber
              value={values[field.key] as number}
              onValueChange={(e) => onChange(field.key, e.value)}
              min={field.min}
              max={field.max}
              step={field.step ?? 1}
              className="w-full"
            />
          )}
          {field.type === "select" && (
            <Dropdown
              value={values[field.key]}
              options={(field.options || VICTOR_TYPES).map((o) => ({ label: o, value: o }))}
              onChange={(e) => onChange(field.key, e.value)}
              className="w-full"
            />
          )}
          {field.type === "checkbox" && (
            <div className="flex items-center gap-2 h-10">
              <Checkbox
                checked={values[field.key] as boolean}
                onChange={(e) => onChange(field.key, e.checked)}
              />
              <span className="text-xs text-slate-500">
                {values[field.key] ? "Enabled" : "Disabled"}
              </span>
            </div>
          )}
        </div>
      ))}
    </div>
  );
}
