declare module "next/font/google" {
  import type { NextFontWithVariable } from "next/dist/compiled/@next/font/dist/types";
  type FontModule = (options: {
    variable?: string;
    subsets?: string[];
    weight?: string | string[];
    style?: string | string[];
    display?: string;
  }) => NextFontWithVariable;
  export const Geist: FontModule;
  export const Geist_Mono: FontModule;
}

declare module "smiles-drawer" {
  const SmilesDrawer: {
    SmiDrawer: new (options?: Record<string, unknown>) => {
      draw(smiles: string, target: Element | string | null, theme?: string, successCb?: (() => void) | null, errorCb?: (() => void) | null): void;
    };
    Parser: {
      parse(smiles: string): unknown;
    };
  };
  export default SmilesDrawer;
}

declare module "next/navigation" {
  export function useParams<T = Record<string, string>>(): T;
  export function useRouter(): {
    push: (url: string) => void;
    replace: (url: string) => void;
    refresh: () => void;
    back: () => void;
    forward: () => void;
    prefetch: (url: string) => void;
  };
  export function usePathname(): string;
  export function useSearchParams(): URLSearchParams;
  export function redirect(url: string): never;
  export function notFound(): never;
}
