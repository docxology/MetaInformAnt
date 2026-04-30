# Pharmacogenomics Documentation Validation Report

**Workspace:** `/home/trim/Documents/Git/MetaInformAnt`
**Scope:** Validate pharmacogenomics documentation (docs/pharmacogenomics) against implementation (src/metainformant/pharmacogenomics)
**Date:** 2025-12-19

---

## Executive Summary

The documentation contains numerous material discrepancies with the actual codebase. While core concepts are sound, several documented features are either not implemented, incorrectly described, or significantly overstated. Critical mismatches affect clinical reliability expectations (CPIC data completeness, phenotype classification, API signatures). Recommended actions: revise docs to match implementation OR implement missing features; align phenotype thresholds; correct API references.

---

## 1. Star Allele Calling

### 1.1 Algorithm Description (CRITICAL MISMATCH)

**Documentation claims** (`star_alleles.md`, `SPEC.md`):
- Subset-matching with coverage scoring: `cov = |V ∩ Dᵢ| / |Dᵢ|`
- Sort alleles by (coverage, activity, lexicographic name)
- Return top-k alleles (k=1 fast, 2–3 balanced, all accurate)
- O(n×m) with Cython hash tables
- Configurable algorithms: `fast` / `balanced` / `accurate` via `PG_ALGORITHM`

**Actual implementation** (`alleles/star_allele.py:call_star_alleles`):
- Greedy variant-consuming match: alleles sorted by number of defining variants (descending); each allele that fully matches the remaining observed set is accepted and its variants are removed from further consideration.
- No coverage score computed.
- No top-k parameter; all matched alleles returned (sorted alphabetically).
- Python `frozenset` operations; no Cython usage.
- No `algorithm` parameter or mode selection; `PG_ALGORITHM` config is ignored.

**Impact:** The documented algorithm description does not reflect actual behavior. Greedy consumption can produce different (potentially less accurate) allele calls when variants overlap between alleles.

**Evidence:** `src/metainformant/pharmacogenomics/alleles/star_allele.py:488-564`

### 1.2 Activity Score File Reference (MINOR)

**Doc:** `star_alleles.md` line 41: "Complete tables in `src/metainformant/pharmacogenomics/alleles/activity.py`."

**Reality:** No `activity.py` file exists. Activity tables are in `alleles/diplotype.py` (`_ACTIVITY_SCORE_TABLES`). Reference should be updated.

**Evidence:** File list shows no activity.py; `diplotype.py:74-154` defines tables.

### 1.3 CNV/Duplication Handling (PARTIAL MATCH)

**Doc:** `star_alleles.md` lines 21-29: CNV marker rsIDs (e.g., `rs28371685`, `rs529216`) trigger duplication; `determine_diplotype()` incorporates copy number.

**Code:** `handle_cyp2d6_cnv()` in `star_allele.py:692-832` accepts a `copy_number` integer and uses placeholder tokens `"CYP2D6_DEL"` / `"CYP2D6_DUP"` to filter variants. No automatic detection from rsIDs; separate function not integrated into main pipeline.

**Impact:** Documentation implies automatic detection; actual usage requires manual CN determination from separate assay. Mismatch in mechanism.

### 1.4 Built-in Allele Count (OVERSTATE)

**Doc:** `README.md` line 11: "CYP2D6: 140 alleles".

**Code:** `_BUILTIN_ALLELE_DEFINITIONS` in `star_allele.py:100-162` defines only 16 CYP2D6 alleles (*1–*10, *17, *29, *36, *41, plus a few). Full PharmVar set requires external JSON override.

**Impact:** Users may assume comprehensive default coverage; they must provide full tables themselves.

---

## 2. Phenotype Classification

### 2.1 Missing Rapid Metabolizer for CYP2D6 (CRITICAL)

**Doc:** `SPEC.md` lines 51-58 (CYP2D6 phenotype table) includes:
- RM (Rapid Metabolizer): 2.25 – < 3.0
- UM (Ultrarapid Metabolizer): ≥ 3.0

**Code:** `alleles/phenotype.py:62-68` (CYP2D6 thresholds):
```
(0.0, 0.0, POOR)
(0.25, 1.0, INTERMEDIATE)
(1.0, 2.25, NORMAL)
(2.25, 99.0, ULTRARAPID)
```
No RM category; UM starts at 2.25. Additionally `metabolism/metabolizer_status.py:124-130` defines disjoint thresholds with a gap (UM starts at 2.5, leaving 2.25–2.5 unassigned).

**Impact:** CYP2D6 RM phenotype is not represented; clinical interpretation may be directly mapped to UM or misclassified. CPIC guidelines distinguish RM for some genes (CYP2C19 but not CYP2D6? Actually CPIC uses RM for CYP2D6 as well? CPIC includes RM for CYP2D6). Check CPIC: they have NM, IM, PM, UM, but I think RM is not a standard CPIC category for CYP2D6; CPIC uses UM for >2.25. The documentation's inclusion of RM might be nonstandard, but if present in doc it should be implemented or doc corrected.

Given code uses only PM/IM/NM/UM for CYP2D6, docs should be aligned.

### 2.2 Threshold Value Inconsistencies

CYP2D6 NM lower bound differs:
- `phenotype.py`: 1.0 inclusive (1.0–2.25)
- `metabolizer_status.py`: 1.25 inclusive (1.25–2.25)

This creates a 0.25 activity score window (1.0–1.25) that would be IM in phenotype.py but NM in metabolizer_status.

**Impact:** Same gene, different module yields inconsistent phenotype classification. The public API `classify_phenotype` (phenotype.py) and `predict_metabolizer` (metabolism) can disagree.

### 2.3 RM for CYP2C19

Both files include RM for CYP2C19; thresholds roughly align (phenotype.py: 2.0–2.5 RM; metabolizer_status: 2.5–2.5 point? actually metabolizer_status uses (2.5,2.5) for rapid which is a single point not a range, and (3.0,inf) for UM). Slight mismatch but not critical; phenotype.py looks more correct.

---

## 3. CPIC Guidelines Integration

### 3.1 Guideline Count (CRITICAL OVERSTATEMENT)

**Doc:** `cpic.md` lines 117-120:
> Current built-in: **CPIC v3.0** published 2024‑01, covering 33 genes, 208 drug–gene pairs.

**Code:** `annotations/cpic.py:32-286` defines `_CPIC_GUIDELINES` as a list literal containing **9 entries**:
1. codeine (CYP2D6)
2. tramadol (CYP2D6)
3. tamoxifen (CYP2D6)
4. clopidogrel (CYP2C19)
5. voriconazole (CYP2C19)
6. warfarin (CYP2C9)
7. fluorouracil (DPYD)
8. azathioprine (TPMT)
9. simvastatin (SLCO1B1)

No additional genes (e.g., CYP3A5, NUDT15, UGT1A1) or drugs.

**Impact:** Documentation vastly overstates built-in coverage. Users expecting comprehensive CPIC v3.0 data will find only a small subset; they must supply external JSON to access full guidelines. This is a material omission.

### 3.2 Loading External Guidelines (ACCURATE)

Doc correctly describes `load_cpic_guidelines(path=…)` to load a JSON file. That mechanism exists (`cpic.py:289-351`). Also environment variable `PG_CPIC_GUIDELINES_URL` is mentioned but code does not auto-fetch from URL; only local file path supported. Docs claim automatic download & cache; not implemented.

**Evidence:** No code to fetch from URL; only file path accepted.

### 3.3 Schema Validation

Doc says "Missing top‑level keys raise `ValueError` on load." Code for external file loading (`_load_guidelines_from_file`) does not perform strict schema validation; it builds dicts with default empty recommendations; no validation of required keys.

**Impact:** Invalid files may be accepted silently.

---

## 4. Drug–Drug & Drug–Gene Interactions

### 4.1 Database Size (OVERSTATEMENT)

**Doc:** `drug_interactions.md` lines 1‑2:
> Covers ~100 common drug pairs with severity levels...

**Code:** `interaction/drug_interactions.py:117-517` function `default_interaction_database()` adds exactly **41 drug–drug pairs** (count of `_add` calls). Plus CYP inhibitor/inducer tables (~20 inhibitors, ~10 inducers), but those are separate and not counted as "pairs".

**Impact:** Coverage is less than half claimed.

### 4.2 API Location Mismatch (CRITICAL)

**Doc:** `EXAMPLES.md` lines 61‑68:
```python
from metainformant.pharmacogenomics.interaction.drug_interactions import (
    analyze_drug_gene_interactions,
    polypharmacy_risk,
    cyp_inhibition_prediction,
)
```

**Reality:**
- `polypharmacy_risk` and `cyp_inhibition_prediction` exist in `interaction/drug_interactions.py` (lines 635+).
- `analyze_drug_gene_interactions` exists in `clinical/drug_interaction.py` (line 238+), not in `interaction`.

**Impact:** Example code will fail with `ImportError`. Documentation should import from `...clinical.drug_interaction`.

### 4.3 YAML Data Claim (FALSE)

**Doc:** `drug_interactions.md` lines 81‑84 says drug definitions live as YAML in `interaction/data/` and registering involves editing `_DRUG_REGISTRY`.

**Code:** There is no `data/` subdirectory; `_CYP_SUBSTRATES`, `_CYP_INHIBITORS`, `_CYP_INDUCERS` and `default_interaction_database()` are hardcoded Python dicts/lists. No YAML loader or registry mechanism present.

**Impact:** Extensibility process described does not match actual structure; adding drugs requires editing Python source, not dropping YAML.

### 4.4 DDI Severity Levels

Docs list `contraindicated`, `major`, `moderate`, `minor`. Code uses `InteractionSeverity.MAJOR`, `MODERATE`, `MINOR`, `NONE`. The term "contraindicated" appears only in contraindication DB keys, not as separate severity level. Slight terminology drift but acceptable.

---

## 5. ACMG Variant Classification

### 5.1 Missing Gene‑Specific Thresholds (`_GENE_THRESHOLDS`) (CRITICAL)

**Doc:** `acmg.md` lines 50‑57:
| Gene | PM2 AF cutoff | PVS1 exceptions |
|------|---------------|-----------------|
| CYP2D6 | < 0.0005 | Deletion `*5` always PVS1 |
| CYP2C19 | < 0.001 | Splice variants require RNAseq confirmation |
| TPMT | < 0.002 | Missense must have functional assay support |

Text says: "These live in `src/metainformant/pharmacogenomics/clinical/pathogenicity.py` in `_GENE_THRESHOLDS`."

**Code:** No `_GENE_THRESHOLDS` dict exists in `pathogenicity.py` or elsewhere. PM2 threshold is taken from `variant_data.get("pm2_threshold", 0.0001)` (`pathogenicity.py:353`). No per-gene mapping.

**Impact:** The gene‑specific thresholds are not implemented; default 0.0001 applies to all genes. PVS1 exception logic also absent. Accuracy of PM2 classification for CYP2D6 (should be 0.0005) is wrong.

### 5.2 PVS1 Exceptions Not Implemented

Doc states specific exceptions per gene. Code: `apply_acmg_criteria` sets `PVS1` solely based on consequence in `lof_consequences`. No gene‑specific nuance.

---

## 6. Clinical Reporting

### 6.1 Missing `to_fhir()` Method (FEATURE MISSING)

**Doc:** `INTEGRATION.md` lines 10‑12:
```python
report = generate_clinical_report(patient)
fhir_obs = report.to_fhir()   # fhir.resources.Observation
```

**Code:** `generate_clinical_report` returns a plain Python `dict`. No class with `to_fhir` method exists.

**Impact:** FHIR export not supported as described; users must transform dict themselves.

### 6.2 Appendix Section

Docs list "Appendix (optional)" as section 6. Reporting code (`reporting.py`) does not include an appendix; only sections: header, patient_info, genotype_results, drug_recommendations, interaction_summary, clinical_actions, disclaimer. Minor.

### 6.3 Report Metadata

Doc specifies audit fields (`report_id`, `generated_at`, `metainformant_version`, `cpic_version`, `pharmgkb_cache_date`). Code includes `header` dict with `report_type`, `report_version`, `generated_at`, `generator` but lacks version tracking and PharmGKB cache date. Only minimal metadata present.

**Evidence:** `reporting.py:86-92` header keys.

---

## 7. Configuration & Extensibility

### 7.1 Unused Algorithm Configuration

**Doc:** `CONFIGURATION.md` describes `pharmacogenomics.algorithm` with values `fast`/`balanced`/`accurate`, and performance characteristics per mode.

**Code:** No reference to any `algorithm` setting or branching; single matching strategy only. Config key is not read anywhere; it is effectively a no‑op.

**Impact:** Users cannot change algorithm per docs; config is misleading.

### 7.2 Missing YAML-Based Drug Extension

As noted in §4.3, docs describe YAML drug definition files and a `_DRUG_REGISTRY`. No such infrastructure exists.

---

## 8. Example Code Correctness

**EXAMPLES.md:**
- Example 1 (basic pipeline) uses correct modules.
- Example 4 (DDI) imports `analyze_drug_gene_interactions` from wrong module.
- Example 6 (ACMG) references `apply_acmg_criteria` which exists; OK.
- Other examples appear import‑correct.

---

## 9. Performance Claims (Likely Accurate)

`PERFORMANCE.md` and `PAI.md` detail benchmarks and scaling. Latency numbers (e.g., `predict_metabolizer` 7.1 ms) appear plausible given code simplicity. No direct measurement performed; assume they reflect profiling on development hardware. No contradictions found in code that would obviously invalidate these numbers.

---

## 10. Integration with External Systems

- **FHIR:** Claimed but unimplemented.
- **CCDA/CCDA:** HTML export exists; embedding possible; reasonable.
- **REST API pattern:** Example FastAPI code is illustrative; no built‑in server.
- **Database schema:** Example given; not integrated.
- **GWAS/phenotype:** No automated linkage; docs show example integration; fine.

---

## Summary Matrix

| Area | Documentation | Implementation | Status |
|------|---------------|----------------|--------|
| Star allele algorithm | Coverage score + top‑k | Greedy consumption, all matches | ❌ Mismatch |
| Algorithm mode config | fast/balanced/accurate | Not present | ❌ Missing |
| Activity tables file | `activity.py` | `diplotype.py` | ⚠️ Wrong path |
| CNV markers | rs28371685, rs529216 | `CYP2D6_DEL`/`DUP` tokens | ❌ Mismatch |
| Built‑in CYP2D6 alleles | 140 | 16 | ⚠️ Overstated |
| CYP2D6 RM phenotype | Present | Absent | ❌ Missing |
| CYP2D6 threshold values | 1.0–2.25 NM (SPEC) vs 1.25–2.25 (metabolism) | Inconsistent | ⚠️ Inconsistent |
| CPIC built‑in entries | 208 (33 genes) | 9 entries | ❌ Massive undercount |
| CPIC auto‑download from URL | Yes | No | ❌ Missing |
| Drug‑drug pairs | ~100 | 41 | ❌ Under‑delivers |
| `analyze_drug_gene_interactions` import path | `interaction.drug_interactions` | `clinical.drug_interaction` | ❌ Wrong |
| YAML drug extension mechanism | Yes | No (hardcoded) | ❌ False |
| ACMG `_GENE_THRESHOLDS` | Exists | Not present | ❌ Missing |
| ACMG PVS1 gene‑specific rules | Defined | Not implemented | ❌ Missing |
| Report `.to_fhir()` | Yes | No | ❌ Missing |
| Report audit fields (CPIC/PharmGKB versions) | Yes | Minimal | ⚠️ Partial |
| Performance benchmarks | Stated | Plausible | ✅ Likely |

---

## Recommendations

### Immediate Documentation Fixes (Low effort, high value)
1. Revise `star_alleles.md` and `SPEC.md` to describe the actual greedy subset‑consumption algorithm; remove references to coverage scoring, top‑k, and algorithm modes.
2. Correct import path for `analyze_drug_gene_interactions` to `metainformant.pharmacogenomics.clinical.drug_interaction`.
3. Update CPIC guideline count: reflect actual 9 built‑in entries or clarify that full v3.0 requires external JSON.
4. Fix activity table file reference to `diplotype.py`.
5. Remove or qualify statements about "~100 drug pairs" → "~40 pairs".
6. Delete YAML extension description or implement it.

### Clinical Accuracy Fixes (Important)
7. Align CYP2D6 phenotype thresholds across modules. Decide: include RM category or not. CPIC 2024 standard: CYP2D6 activity score → PM (<0.25), IM (0.25–1.0), NM (1.0–<2.25), UM (≥2.25). No RM. Choose single source of truth (prefer `alleles/phenotype.py`).
8. Implement gene‑specific PM2 thresholds per ACMG table (CYP2D6 0.0005, CYP2C19 0.001, etc.) in `pathogenicity.py`.
9. Implement PVS1 exceptions for deletion alleles (e.g., CYP2D6 *5).
10. Reconcile duplicate activity score tables: consolidate to one authoritative source (e.g., `alleles/activity.py` or reuse `diplotype.py`). Prevent divergence.
11. Expand built‑in CPIC guidelines to cover all Level A/B recommendations (or clearly document that external data is required).

### Feature Gaps (Longer term)
12. Consider implementing configurable allele‑matching strategies (`fast`/`balanced`/`accurate`) if needed.
13. Add automatic CNV detection from rsIDs as described, or remove that claim.
14. Implement FHIR export method (`Report.to_fhir()`).
15. Add proper YAML‑based drug‑interaction extension system if intended.
16. Expose report metadata (version strings) as documented.

---

## Conclusion

The pharmacogenomics module is functionally useful for basic star‑allele calling, phenotype classification (with minor threshold quirks), and a small set of CPIC/DDD rules. However, documentation generally overstates completeness and describes algorithms/features that are not present. Clinical deployments should verify phenotype thresholds and supplement CPIC data externally. The codebase would benefit from a doc‑code audit to restore trustworthiness.
