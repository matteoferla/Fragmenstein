# Plan: Expand Analog Search Sources (Similar Search Space)

## Context

The Similars step currently only supports SmallWorld (sw.docking.org) as the search backend. This is limiting because:
- SmallWorld has SSL/connectivity issues (broken TLS 1.3 on their server)
- Scientists often already have their own compound lists from vendor catalogs, literature, or internal libraries
- Different databases cover different chemical space (purchasable vs bioactive vs patent-derived)
- SmallWorld is graph-edit-distance based — fingerprint similarity searches (PubChem, ChEMBL) find different types of analogs

The Similars page should offer **multiple modes** for providing analogs to the Place step.

---

## Proposed Modes

### Mode 1: SmallWorld Search (existing)
Already implemented. Searches sw.docking.org for graph-edit-distance neighbors in commercial databases (Enamine REAL, ZINC, Mcule, etc.).

**No changes needed** — already works.

---

### Mode 2: Manual SMILES Input (paste)
User pastes SMILES directly into a text area. Format: one compound per line, `SMILES name` (tab or space separated). Name is optional — auto-generate if missing.

**Why:** Fastest way for a chemist to test a handful of specific compounds. No file needed.

**Backend:**
- New endpoint: `POST /api/sessions/{session_id}/similars/manual`
- Request body: `{ smiles_text: string, combine_job_id?: string }`
- Service: Parse text, validate each SMILES via RDKit, build DataFrame with `smiles`, `name` columns, save as similars results (same format as SmallWorld output, just fewer columns)
- The place step's `build_place_queries_from_similars` already handles DataFrames with just `smiles` + `name` — it falls back to all session hits when `hit_mols` column is absent

**Frontend:**
- Add a TabView on the similars page: "SmallWorld" | "Paste SMILES" | "Upload File"
- Paste tab: `<InputTextarea>` with placeholder showing format example, a "Load" button
- No background job needed — parsing is instant, just POST and get results back synchronously
- After loading, show the same `SimilarsTable` with 2D depictions

**Files to modify:**
- `web/backend/app/routers/similars.py` — new endpoint
- `web/backend/app/services/similars_service.py` — new `parse_manual_smiles()` function
- `web/frontend/src/app/sessions/[id]/similars/page.tsx` — add TabView with paste mode

---

### Mode 3: CSV/Excel Upload
User uploads a CSV, TSV, or Excel (.xlsx) file containing SMILES and compound identifiers.

**Why:** Scientists export compound lists from vendor catalogs (Enamine, ChemSpace, eMolecules), internal databases, or spreadsheets. This is the most common workflow for experienced users.

**Expected file formats:**
- CSV/TSV with header row. Must have a column named `smiles` or `SMILES`. Optional columns: `name`, `id`, `vendor_id`, or similar for the compound identifier. Any extra columns (MW, price, vendor) are preserved and shown in the results table.
- Excel (.xlsx) — same column expectations, reads first sheet.

**Backend:**
- New endpoint: `POST /api/sessions/{session_id}/similars/upload` (multipart file upload)
- Service: Read CSV/TSV via `pd.read_csv()` or Excel via `pd.read_excel()`. Auto-detect the SMILES column (look for column named `smiles`, `SMILES`, `canonical_smiles`, `Smiles`). Validate each SMILES with RDKit. Auto-detect name column (`name`, `Name`, `id`, `ID`, `compound_id`, `vendor_id`). Save as similars results file.
- Need `openpyxl` dependency for Excel support.

**Frontend:**
- Upload tab in the TabView: `<FileUpload>` component accepting `.csv, .tsv, .xlsx`
- After upload, show results in `SimilarsTable`
- Show a summary: "Loaded 150 compounds from catalog.xlsx (3 invalid SMILES skipped)"

**Files to modify:**
- `web/backend/app/routers/similars.py` — new upload endpoint
- `web/backend/app/services/similars_service.py` — new `parse_uploaded_file()` function
- `web/frontend/src/app/sessions/[id]/similars/page.tsx` — add upload tab
- `web/backend/requirements.txt` — add `openpyxl`

---

### Mode 4: PubChem Similarity Search (future)
Search PubChem's 100M+ compound database using fingerprint similarity (Tanimoto). Free, no auth required, reliable servers.

**Why:** Different chemical space from SmallWorld. PubChem includes bioactive compounds, natural products, and approved drugs — useful for scaffold hopping and repurposing. No SSL issues (Google infrastructure).

**API:** PubChem PUG REST — `GET /rest/pug/compound/fastsimilarity_2d/smiles/{smiles}/property/CanonicalSMILES,IUPACName,MolecularWeight/JSON?Threshold=80`

**Backend:**
- New service: `pubchem_service.py` — search PubChem for each merger SMILES, collect results
- Uses `requests.get()` (no SSL issues — PubChem uses standard TLS)
- Returns DataFrame with `smiles`, `name` (CID), `tanimoto`, `molecular_weight`, `iupac_name`

**Frontend:**
- Add "PubChem" tab to the TabView
- Config: Tanimoto threshold slider (70-100%), max results per query
- Results in same `SimilarsTable` format

**Files to create:**
- `web/backend/app/services/pubchem_service.py` — new
- Extend `web/backend/app/routers/similars.py` — new endpoint
- Extend frontend similars page — new tab

---

## Implementation Priority

| Mode | Effort | Value | Priority |
|------|--------|-------|----------|
| Manual SMILES paste | S | High — instant, no external deps | 1st |
| CSV/Excel upload | S-M | High — most common real workflow | 2nd |
| PubChem similarity | M | Medium — free alternative search | 3rd |
| SmallWorld | Done | Done | Done |

---

## Shared Design Principle

All modes produce the **same output format**: a DataFrame saved to disk with at minimum `smiles` and `name` columns. This is stored as the similars job result and consumed by `build_place_queries_from_similars()` in the Place step. The Place step doesn't know or care which source the analogs came from.

The frontend `SimilarsTable` component already handles varying columns gracefully — SmallWorld-specific columns (`topodist`, `ecfp4`, etc.) only appear when present.

---

## Frontend UX: Tab-Based Mode Selector

```
+--────────────+──────────────+──────────────+──────────────+
|  SmallWorld  | Paste SMILES | Upload File  |   PubChem    |
+──────────────+──────────────+──────────────+──────────────+
|                                                           |
|  [Mode-specific config/input area]                        |
|                                                           |
|  [Action button: Search / Load / Upload]                  |
|                                                           |
+───────────────────────────────────────────────────────────+
|                                                           |
|  [SimilarsTable -- shared across all modes]               |
|  [2D structures, SMILES, vendor ID, metrics]              |
|                                                           |
|  [Detail panel on row select]                             |
|                                                           |
+───────────────────────────────────────────────────────────+
```

The results table and detail panel below are identical regardless of source mode. Only the input area above the table changes per tab.

---

## Backend Endpoints Summary

| Method | Path | Mode |
|--------|------|------|
| `POST` | `/api/sessions/{id}/similars` | SmallWorld (existing) |
| `POST` | `/api/sessions/{id}/similars/manual` | Paste SMILES |
| `POST` | `/api/sessions/{id}/similars/upload` | CSV/Excel upload |
| `POST` | `/api/sessions/{id}/similars/pubchem` | PubChem search (future) |

All endpoints save results to the same path and create a job record so the Place step can reference it via `source_job_id`.

---

## Verification Plan

1. **Paste mode**: Paste 5 SMILES, verify they appear in SimilarsTable with 2D depictions, proceed to Place
2. **Upload mode**: Upload a CSV with 20 compounds (SMILES + name columns), verify parsing, invalid SMILES reported and skipped
3. **Upload Excel**: Same with .xlsx file
4. **Place integration**: After loading analogs from any mode, click "Place Analogs", verify the Place step picks them up correctly
5. **Mixed workflow**: Run SmallWorld, then switch to paste mode and add more compounds, verify both sets are available
