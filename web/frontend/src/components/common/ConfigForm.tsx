"use client";

import { InputText } from "primereact/inputtext";
import { InputNumber } from "primereact/inputnumber";
import { Dropdown } from "primereact/dropdown";
import { Checkbox } from "primereact/checkbox";
import { VICTOR_TYPES } from "@/lib/constants";

export interface ConfigField {
  key: string;
  label: string;
  type: "text" | "number" | "select" | "checkbox";
  options?: string[];
  optionDescs?: Record<string, string>;
  description?: string;
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
      {fields.map((field) => {
        const selectedValue = values[field.key];
        const optionDesc = field.optionDescs && typeof selectedValue === "string"
          ? field.optionDescs[selectedValue]
          : null;

        return (
          <div key={field.key} className="flex flex-col gap-1.5">
            <label className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">
              {field.label}
            </label>
            {field.type === "text" && (
              <InputText
                value={(selectedValue as string) || ""}
                onChange={(e) => onChange(field.key, e.target.value)}
                className="w-full"
              />
            )}
            {field.type === "number" && (
              <InputNumber
                value={selectedValue as number}
                onValueChange={(e) => onChange(field.key, e.value)}
                min={field.min}
                max={field.max}
                step={field.step ?? 1}
                className="w-full"
              />
            )}
            {field.type === "select" && (
              <Dropdown
                value={selectedValue}
                options={(field.options || VICTOR_TYPES).map((o) => ({ label: o, value: o }))}
                onChange={(e) => onChange(field.key, e.value)}
                className="w-full"
              />
            )}
            {field.type === "checkbox" && (
              <div className="flex items-center gap-2 h-10">
                <Checkbox
                  checked={selectedValue as boolean}
                  onChange={(e) => onChange(field.key, e.checked)}
                />
                <span className="text-xs text-slate-500">
                  {selectedValue ? "Enabled" : "Disabled"}
                </span>
              </div>
            )}
            {/* Description: static or dynamic based on selected option */}
            {optionDesc && (
              <p className="text-[11px] text-slate-400 leading-relaxed -mt-0.5">{optionDesc}</p>
            )}
            {field.description && !optionDesc && (
              <p className="text-[11px] text-slate-400 leading-relaxed -mt-0.5">{field.description}</p>
            )}
          </div>
        );
      })}
    </div>
  );
}
