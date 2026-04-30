# GWAS Pipeline Documentation Validation Report

**Date**: 2025-10-30
**Scope**: Validate all QC, association, and visualization methods in docs/gwas/ against src/metainformant/gwas/ implementations
**Status**: ⚠️ Multiple discrepancies found (see prioritized fixes below)

---

## Executive Summary

This review compared the GWAS pipeline documentation (`docs/gwas/`) with the actual implementation (`src/metainformant/gwas/`). The core pipeline is largely functional with robust QC, association testing, and visualization. However, **several critical gaps** exist between documentation and implementation, including:

- **CLI command documented but not implemented** (`python -m metainformant gwas run`)
- **ADMIXTURE computational method missing** (only visualization of pre-computed results exists)
- **QC parameter default mismatch** (max_missing: documented 0.05, implemented 0.1)
- **QC parameter name inconsistency** (`hwe_pval` vs `min_hwe_p`)
- **min_call_rate filter not implemented** (mapped in config but never applied)
- **Permutation testing claimed complete but not implemented**

---

## ✅ VERIFIED: Working Components

### 1. VCF Processing and QC Pipeline
**Status**: ✅ Implemented and matches core functionality

**Implementation**: `src/metainformant/gwas/analysis/quality.py`
- `parse_vcf_full()`: Full VCF parsing with NumPy optimization
- `filter_by_maf()`: MAF filtering
- `filter_by_missing()`: Missing data filtering
- `test_hwe()`: Hardy-Weinberg equilibrium chi-square test
- `apply_qc_filters()`: Comprehensive QC pipeline

**Documentation match**: Core QC steps (MAF, missingness, HWE) are correctly documented in `workflow.md` and `config.md`.

---

### 2. Association Testing
**Status**: ✅ Implemented

**Implementation**: `src/metainformant/gwas/analysis/association.py`
- `association_test_linear()`: Linear regression with covariate adjustment
- `run_linear_model_gwas()`: Batch linear GWAS (vectorized)
- `association_test_logistic()`: Logistic regression with odds ratio
- `run_logistic_model_gwas()`: Batch logistic GWAS

**Mixed Model**: ✅ Implemented but **documentation outdated**
- `src/metainformant/gwas/analysis/mixed_model.py`: Full EMMA algorithm (Kang et al. 2008)
- Functions: `association_test_mixed()`, `run_mixed_model_gwas()`
- Eigendecomposition, REML variance component estimation, SNP testing
- Documentation `review.md` incorrectly claims MLM is "NOT IMPLEMENTED"

---

### 3. Multiple Testing Correction
**Status**: ✅ Implemented as documented

**Implementation**: `src/metainformant/gwas/analysis/correction.py`
- `bonferroni_correction()`: Family-wise error rate control
- `fdr_correction()`: Benjamini-Hochberg FDR
- `genomic_control()`: Lambda GC calculation

**Note**: Documentation `review.md` claims "Permutation Testing: Complete" but permutation testing is **NOT implemented** anywhere in the codebase.

---

### 4. Visualization
**Status**: ✅ Implemented

**Implementation**: `src/metainformant/gwas/visualization/general.py`
- `manhattan_plot()`: Genome-wide association plot with chromosome coloring, significance thresholds, gene annotations
- `qq_plot()`: Enhanced Q-Q plot with λ_GC, 95% CI, GC-corrected overlay, KS test
- `generate_all_plots()`: Batch plot generation

**Additional visualizations** in `visualization/population/`:
- PCA scatter plots
- Admixture plots (visualization only – reads pre-computed admixture proportions)
- Geographic distribution

---

### 5. Population Structure (PCA and Kinship)
**Status**: ✅ PCA and kinship implemented; ADMIXTURE computation missing

**Implementation**: `src/metainformant/gwas/analysis/structure.py`
- `compute_pca()`: SVD-based PCA with missing data imputation
- `compute_kinship_matrix()`: Supports VanRaden, Astle, Yang, and IBS methods
- `estimate_population_structure()`: Orchestrates PCA + kinship

**Gap**: Documentation mentions "ADMIXTURE" as a population structure method, but the code only includes **visualization of admixture proportions** (`visualization/population/population_admixture.py`). There is no function to actually run ADMIXTURE software to compute ancestry proportions. The documentation should clarify that ADMIXTURE results must be computed externally.

---

### 6. Fine-Mapping and Colocalization
**Status**: ✅ Implemented

**Implementation**:
- `src/metainformant/gwas/finemapping/credible_sets.py`: SuSiE regression, credible sets, Bayes factors
- `src/metainformant/gwas/finemapping/colocalization.py`: eQTL colocalization ( coloc, CLPP, regional)

---

### 7. Heritability Estimation
**Status**: ✅ Implemented

**Implementation**: `src/metainformant/gwas/heritability/estimation.py`
- LD Score Regression (`estimate_h2_ldsc()`)
- Partitioned heritability
- Genetic correlation
- GREML and Haseman-Elston regression

---

### 8. Benchmarking Module
**Status**: ✅ Implemented and matches documentation

**Implementation**: `src/metainformant/gwas/analysis/benchmarking.py`
- `benchmark_subset_run()`: Pilot timing
- `extrapolate_full_genome_time()`: Scaling extrapolation
- `scaling_model()`: Complexity models (O(n·m), O(n²·m), O(m·k²), etc.)

**Complexity models in code match** `docs/gwas/BENCHMARKING.md` Table 1 exactly.

---

## ❌ CRITICAL DISCREPANCIES

### Issue 1 — CLI Command Not Implemented
**Severity**: 🔴 Critical (breaks user workflows)

**Documentation claim** (throughout `docs/gwas/`):
```bash
python -m metainformant gwas run --config config/gwas/gwas_template.yaml
```

**Actual implementation**: `src/metainformant/__main__.py` only supports `metainformant gwas info`. The `run` subcommand is **not implemented**. The GWAS handler returns 1 for any command other than `info`.

**Workaround**: Actual execution uses standalone scripts in `scripts/gwas/`:
```bash
uv run python scripts/gwas/run_pbarbatus_gwas.py
uv run python scripts/gwas/run_amellifera_gwas.py
```

**Fix priority**: HIGH – Documentation must be corrected or CLI must be implemented.

---

### Issue 2 — ADMIXTURE Misrepresented as Computational Method
**Severity**: 🟠 High (misleads users)

**Documentation claim** (`docs/gwas/workflow.md`, `AGENTS.md`):
> Population structure: PCA, ADMIXTURE

**Implementation**: Only PCA is computed. Admixture visualization (`visualization/population/population_admixture.py`) only **reads** pre-computed admixture proportions from a file (e.g., from the ADMIXTURE software). There is **no function** to run ADMIXTURE or compute ancestry proportions.

**Impact**: Users cannot run ADMIXTURE through the pipeline; they must run external software manually.

**Fix priority**: HIGH – Either implement ADMIXTURE integration or clarify documentation.

---

### Issue 3 — QC Default Mismatch: `max_missing`
**Severity**: 🟡 Medium (affects results reproducibility)

| Parameter | Documentation (`config.md`) | Implementation (`quality.py:417`) |
|-----------|----------------------------|-----------------------------------|
| max_missing | 0.05 | **0.1** |

The default maximum missing genotype proportion differs. The template in `src/metainformant/gwas/data/config.py` uses 0.05, but the actual `apply_qc_filters()` defaults to 0.1 when not provided.

**Risk**: Users getting different results depending on whether they rely on config template vs code defaults.

**Fix priority**: MEDIUM – Align code default with documented/template default (0.05) or document the actual default.

---

### Issue 4 — QC Parameter Name Inconsistency
**Severity**: 🟡 Medium (confusing configuration)

**Documentation** (`config.md` line 102): `hwe_pval`
**Implementation** (`quality.py`): `min_hwe_p`

Although `workflow_config.py` correctly maps `hwe_pval` → `min_hwe_p`, direct API users must use `min_hwe_p`. The parameter name in the `apply_qc_filters()` docstring says `min_hwe_p` but the configuration layer hides this.

**Risk**: API users reading code docstrings see different names than config-file users.

**Fix priority**: MEDIUM – Standardize on `min_hwe_p` internally and document both names clearly, or rename to `hwe_pval` everywhere.

---

### Issue 5 — `min_call_rate` Configuration Mapped but Not Applied
**Severity**: 🔴 Critical (filter silently skipped)

**Evidence**:
- `workflow_config.py` lines 148–149: Maps `min_call_rate` from config to QC parameters
- `quality.py` `apply_qc_filters()`: Never reads `min_call_rate` from `qc_config`; only uses `min_maf`, `max_missing`, `min_hwe_p`

**Impact**: Users configuring `min_call_rate` believe it will be enforced, but it is silently ignored. No `filter_by_call_rate()` function exists.

**Fix priority**: HIGH – Either implement call-rate filtering or remove the configuration option.

---

### Issue 6 — Permutation Testing Claimed but Not Implemented
**Severity**: 🟠 High (documentation error)

**Documentation** (`docs/gwas/review.md` line 69):
> **Permutation Testing**: Complete (basic implementation)

**Reality**: No permutation testing code found in `src/metainformant/gwas/analysis/correction.py` or elsewhere.

**Fix priority**: HIGH – Remove claim or implement permutation-based multiple testing.

---

### Issue 7 — Missing `min_qual` Filter Implementation
**Severity**: 🟡 Medium (incomplete QC)

**Documentation** (`config.md` line 104): `min_qual: 30.0` documented as QC filter.

**Implementation**: `apply_qc_filters()` extracts `min_qual` from config (line 400 docstring) but never applies it. No code filters variants by QUAL score.

**Fix priority**: MEDIUM – Implement quality score filtering or remove from documentation.

---

### Issue 8 — Missing `exclude_indels` Filter Implementation
**Severity**: 🟡 Medium (incomplete QC)

Similar to `min_qual`: documented but not applied in `apply_qc_filters()`.

**Fix priority**: MEDIUM – Implement indel exclusion or remove from documentation.

---

## ⚠️ Other Discrepancies

### Discrepancy 9 — `run_gwas()` Function Location
**Documentation** (`README.md` line 56 table): Lists `run_gwas()` as top-level function.
**Implementation**: `run_gwas()` exists in `workflow/workflow_execution.py` but is not imported at package level. Users must import from `metainformant.gwas.workflow` not `metainformant.gwas`.

**Fix**: Update documentation table to show correct import path.

---

### Discrepancy 10 — Visualization Module Names in Docs
**Documentation** (`visualization_gallery.md`) references:
- `visualization_genome.py`
- `visualization_statistical.py`

**Reality**: The module is `visualization/general.py` (not split into genome/statistical files). The gallery incorrectly suggests separate modules.

**Fix**: Correct module names in gallery documentation.

---

### Discrepancy 11 — `generate_all_plots()` Duplicate Implementations
Two `generate_all_plots()` functions exist:
- `visualization/general.py`: line 1209
- `visualization/interactive/suite.py`: line 57

This may be intentional (separate batch vs interactive) but could cause confusion.

---

### Discrepancy 12 — Missing `manhattan_plot()` Parameter in Docs
Documentation shows usage:
```python
manhattan_plot(results="results.tsv", ...)
```

**Implementation** expects `results` as list of dicts or single dict, **not** a file path. The documentation example is wrong – it shows passing a filename string where a data structure is expected.

**Fix**: Correct example to load data first, then pass results dict.

---

### Discrepancy 13 — `qq_plot()` Parameter Mismatch
**Documentation** (`visualization_gallery.md` line 99–104):
```python
qq_plot(results="results.tsv", show_ci=True, show_lambda_gc=True)
```

**Implementation** (`general.py` line 337):
- Accepts `p_values` as list, not `results` filename
- No `show_ci` or `show_lambda_gc` parameters (CI and GC are automatic)

**Fix**: Correct documentation to match actual function signature.

---

### Discrepancy 14 — Workflow Step Names vs Benchmarking Names
`workflow.md` Table lists step: `multiple_testing`
`benchmarking.md` scaling models uses: `multiple_testing_correction`

Name mismatch could cause custom model overrides to fail.

---

### Discrepancy 15 — Association Testing Output Fields
**Documentation** (`workflow.md` line 292):
> Outputs: beta, se, p-value, R² (for linear), odds ratio (for logistic)

**Implementation**: Both linear and logistic return `beta`, `se`, `p_value`, plus specific stats. Correct.

However, `run_linear_model_gwas()` also includes `maf` computed per variant. Not documented but useful.

---

## 📋 CONFIGURATION VALIDATION

### Species-Specific Configs
**Files checked**:
- `docs/gwas/amellifera_config.md` – Apis mellifera
- `docs/gwas/pbarbatus_config.md` – Pogonomyrmex barbatus

**Assessment**: ⚠️ `amellifera_config.md` states "20 principal components" but this is not enforced in code; it's a recommendation, not a default. The default `n_components` is 10 in `compute_pca()`.

**Action**: Document that 20 PCs is recommended for honeybee, not the default.

---

### Template Config vs Implementation
**Template** (`src/metainformant/gwas/data/config.py` `create_gwas_config_template()`):
- maf_threshold: 0.01 ✓
- missing_threshold: 0.05 ✓
- hwe_p_threshold: 1e-6 ✓
- pca_components: 10 ✓

These match the documented template values. However, the code's runtime default in `apply_qc_filters()` for `max_missing` is 0.1, not 0.05. This is a **default value divergence** between the template (what users copy) and the runtime fallback (what happens if they omit the key).

---

## 🎯 PRIORITIZED FIXES

### 🔴 CRITICAL (Must Fix Before Release)

1. **Fix CLI command or documentation**
   - Option A: Implement `metainformant gwas run` subcommand in `__main__.py` that calls `execute_gwas_workflow()`
   - Option B: Remove all mentions of `metainformant gwas run` from docs and direct users to `scripts/gwas/` pipeline scripts
   - **Impact**: Users cannot run GWAS as documented

2. **Implement or remove `min_call_rate` filter**
   - Either add `filter_by_call_rate()` to `quality.py` and apply in `apply_qc_filters()`
   - Or remove `min_call_rate` from `config.md`, `workflow_config.py` mapping, and template config
   - **Impact**: Configuration silently ignored

3. **Implement or remove `min_qual` and `exclude_indels` filters**
   - Add quality filtering and indel exclusion to `apply_qc_filters()`
   - Or remove from documentation and config schema
   - **Impact**: Promised QC not applied

4. **Remove or implement permutation testing claim**
   - Implement permutation-based p-value adjustment if claimed
   - Or remove "Permutation Testing: Complete" from `review.md`
   - **Impact**: Feature falsely advertised

---

### 🟠 HIGH (Fix Before Next Major Update)

5. **Resolve ADMIXTURE documentation ambiguity**
   - Add documentation that ADMIXTURE (software) must be run externally; only visualization of results is provided
   - Consider renaming "ADMIXTURE plots" to "Ancestry proportion plots" to avoid implying software integration
   - **Impact**: Users expect built-in ancestry estimation

6. **Fix `max_missing` default divergence**
   - Change `apply_qc_filters()` line 417 default from 0.1 to 0.05 to match template and documentation
   - Or document that 0.1 is the code default and explain rationale
   - **Impact**: Reproducibility – different default values

7. **Standardize HWE parameter naming**
   - Choose one: `min_hwe_p` (code) or `hwe_pval` (config docs)
   - Update all code, config schema, and documentation consistently
   - **Impact**: API clarity

8. **Correct `review.md` on MLM status**
   - Update to reflect that EMMA mixed model **is implemented** in `mixed_model.py`
   - **Impact**: False negative assessment

---

### 🟡 MEDIUM (Improve Quality)

9. **Fix visualization function signatures in `visualization_gallery.md`**
   - Update examples to pass data structures, not filenames
   - Use correct parameter names (`p_values` not `results`; no `show_ci`)

10. **Document correct import paths**
    - Update `README.md` table to show full import paths (e.g., `from metainformant.gwas.analysis.association import run_linear_model_gwas`)
    - Avoid suggesting top-level imports that don't exist

11. **Standardize workflow step names**
    - Use consistent naming: either `multiple_testing` or `multiple_testing_correction` everywhere
    - Update `workflow.md`, `BENCHMARKING.md`, and code comments

12. **Clarify amellifera PCA recommendation**
    - Mark "20 PCs" as recommended tuning, not default
    - Show config snippet that sets `structure.n_components: 20`

---

## 📊 DETAILED COMPONENT MATRIX

| Component | Implemented | Documented | Match | Notes |
|-----------|-------------|------------|-------|-------|
| VCF parsing | ✅ | ✅ | ✅ | `parse_vcf_full()` matches docs |
| MAF filter | ✅ | ✅ | ✅ | `filter_by_maf()`, default 0.01 |
| Missingness filter | ✅ | ✅ | ⚠️ | Default mismatch: doc 0.05 vs code 0.1 |
| HWE test | ✅ | ✅ | ⚠️ | Param name: `min_hwe_p` vs `hwe_pval` |
| Quality filter | ❌ | ✅ | ❌ | `min_qual` documented but not applied |
| Indel exclusion | ❌ | ✅ | ❌ | Not implemented |
| Call rate filter | ❌ | ✅ | ❌ | `min_call_rate` mapped but unused |
| PCA | ✅ | ✅ | ✅ | SVD-based, default 10 PCs |
| Kinship (VanRaden/Astle/Yang) | ✅ | ✅ | ✅ | All 4 methods (plus IBS) |
| ADMIXTURE computation | ❌ | ✅ | ❌ | Only visualization, no integration |
| Linear regression | ✅ | ✅ | ✅ | Batch and single-variant |
| Logistic regression | ✅ | ✅ | ✅ | With odds ratio |
| Mixed model (EMMA) | ✅ | ❌ | ⚠️ | Implemented but `review.md` says not |
| Bonferroni correction | ✅ | ✅ | ✅ | |
| FDR (BH) | ✅ | ✅ | ✅ | |
| Genomic control | ✅ | ✅ | ✅ | λ_GC calculation |
| Permutation testing | ❌ | ✅ | ❌ | Claimed but absent |
| Manhattan plot | ✅ | ✅ | ⚠️ | API mismatch in docs example |
| Q-Q plot | ✅ | ✅ | ⚠️ | API mismatch in docs example |
| Regional plot | ✅ | ✅ | ✅ | |
| Admixture plot | ✅ (viz only) | ✅ | ⚠️ | Misleading – needs external ADMIXTURE |
| SuSiE fine-mapping | ✅ | ✅ | ✅ | |
| Colocalization (coloc) | ✅ | ✅ | ✅ | |
| LD pruning | ✅ | Not prominently documented | ⚠️ | Implemented, missing from main docs |
| Heritability (LDSC) | ✅ | ✅ | ✅ | |
| Benchmarking module | ✅ | ✅ | ✅ | Complexity models match |
| CLI `gwas run` | ❌ | ✅ | ❌ | Not implemented in `__main__.py` |
| Python API `execute_gwas_workflow()` | ✅ | ✅ | ✅ | |

---

## 📁 FILE-BY-FILE VALIDATION

### Documentation Files
| File | Status | Issues |
|------|--------|--------|
| README.md | ⚠️ Partial | API examples misleading (filename vs data) |
| workflow.md | ⚠️ Partial | CLI command wrong; ADMIXTURE oversold |
| config.md | ⚠️ Partial | Defaults mismatch; unused params listed |
| visualization_gallery.md | ⚠️ Partial | Incorrect module names; wrong function signatures |
| real_data_acquisition.md | ✅ Good | Accurate external tool instructions |
| amellifera_config.md | ⚠️ Partial | PC count presented as default, not config |
| pbarbatus_config.md | ⚠️ Not checked fully | Assumed similar issues |
| BENCHMARKING.md | ✅ Good | Matches implementation |
| assessment.md | Not provided | – |
| verification_report.md | ✅ Good | Reflects test success, but may overclaim |
| review.md | ❌ Outdated | Claims MLM missing, permutation present – both wrong |

---

### Source Code Files
| File | Status | Notes |
|------|--------|-------|
| analysis/quality.py | ⚠️ Partial | Missing min_qual, exclude_indels, min_call_rate |
| analysis/association.py | ✅ Good | Linear, logistic both solid |
| analysis/mixed_model.py | ✅ Good | Full EMMA implementation |
| analysis/correction.py | ✅ Good | Missing permutation (not claimed in code) |
| analysis/benchmarking.py | ✅ Good | Accurate |
| analysis/ld_pruning.py | ✅ Good | Working but not prominently documented |
| analysis/structure.py | ✅ Good | PCA + 4 kinship methods |
| visualization/general.py | ✅ Good | manhattan_plot, qq_plot, generate_all_plots |
| visualization/population/ | ⚠️ Partial | Admixture plot = visualization only |
| workflow/workflow_config.py | ⚠️ Partial | Maps unused params (min_call_rate, min_qual, exclude_indels) |
| workflow/workflow_execution.py | ✅ Good | `run_gwas()` and `execute_gwas_workflow()` exist |
| `__main__.py` | ❌ Incomplete | GWAS only supports `info`, not `run` |

---

## 🔧 RECOMMENDED ACTIONS

### Immediate (Pre-Release)
1. **Choose CLI strategy**: Either implement `gwas run` subcommand in `__main__.py` OR remove all `python -m metainformant gwas run` examples from documentation and point to `scripts/gwas/` instead.
2. **Remove or implement unused QC parameters**: Decide on `min_call_rate`, `min_qual`, `exclude_indels` – implement or prune from config schema and docs.
3. **Fix `max_missing` default**: Align code default (0.1) with template (0.05) OR clearly document the 0.1 default and why.
4. **Standardize parameter names**: Use `min_hwe_p` everywhere; update config docs to match or vice versa.
5. **Update `review.md`**: Correct MLM status (implemented) and remove permutation testing claim.

### Short-term (Next Sprint)
6. Clarify ADMIXTURE's role: document that only plotting pre-computed ancestry is supported, or integrate a simple ancestry estimation method (e.g., NGSAdmix REST wrapper, or basic structure-like algorithm).
7. Add LD pruning to main workflow docs (it exists but is not highlighted).
8. Fix visualization examples to show correct API usage (pass data, not filenames).
9. Unify step naming (`multiple_testing` vs `multiple_testing_correction`) across all docs and code.
10. Document actual import paths in README table.

### Long-term (Future Enhancement)
11. Consider adding a unified `metainformant gwas run` CLI that orchestrates the pipeline (wrap the scripts/ logic).
12. Expand fine-mapping with conditional analysis and credible set refinement.
13. Add permutation-based multiple testing correction.
14. Implement quality filter on variant QUAL scores.
15. Implement indel exclusion filter.

---

## ✅ POSITIVE FINDINGS

- Core QC pipeline robust: MAF, missingness, HWE all correctly implemented and tested.
- Association testing comprehensive: linear, logistic, mixed models (EMMA) all present.
- Multiple testing correction complete: Bonferroni, FDR, genomic control.
- Visualization suite rich: Manhattan, Q-Q, regional, admixture, PCA.
- Fine-mapping and colocalization present and well-structured.
- Heritability estimation comprehensive (LDSC, GREML, HE regression).
- Benchmarking module accurate and matches documented complexity models.
- Species-specific configs (amellifera, pbarbatus) reflect domain knowledge.
- Code quality high: type hints, graceful fallbacks, detailed logging, thorough docstrings.

---

## CONCLUSION

The GWAS pipeline is **largely functional and well-engineered** with strong core methods. However, **documentation is out of sync** with implementation in several areas, and **some promised features are missing** (CLI `gwas run`, QC filters, ADMIXTURE integration). The most critical issues are:

1. **Documented CLI does not exist** – users cannot run `python -m metainformant gwas run` as shown throughout the docs.
2. **ADMIXTURE misrepresented** – only visualization exists, not computation.
3. **Silent configuration failures** – `min_call_rate`, `min_qual`, `exclude_indels` are accepted but ignored.

These issues undermine user trust and reproducibility. Addressing the **🔴 Critical** fixes should be prioritized before promoting the GWAS module to new users.

---

**Report generated by**: Hermes Agent (Nous Research)
**Validation method**: Systematic comparison of documentation files against source code, function signatures, default parameters, and control flow.
