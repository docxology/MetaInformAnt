> Historical snapshot: retained for provenance. Current code, tests, and domain docs are the source of truth.

# MetaInformAnt Documentation Audit — Executive Summary Dashboard

> Historical snapshot: this dashboard is retained for provenance and may not
> describe the current checkout. Regenerate current verification outputs under
> `output/`; the 2026-05-25 stabilization pass confirmed test collection and the
> local non-network/non-external test suite.

**Report Date:** April 29, 2026
**Auditor:** Hermes Agent (Nous Research)
**Scope:** Full documentation audit across 17+ module validation reports
**Dashboard Purpose:** Leadership review — at-a-glance health assessment, critical risks, and actionable roadmap

---

## SNAPSHOT

| Metric | Value | Status |
|--------|-------|--------|
| Composite Documentation Accuracy | 72.8% | ⚠️ Moderate |
| Modules Audited | 14 | — |
| Total Cross-Code Violations | 3,141 | 🔴 Critical |
| Broken Internal Links | 288 | 🔴 Critical |
| Broken External Links | 14 | 🟠 High |
| Missing/Incorrect Examples | 60+ | 🔴 Critical |
| P0 Critical Issues | 12 | 🔴 Release-blockers |
| P1 High Issues | 28 | 🟠 Major |
| P2 Medium Issues | 47 | 🟡 Significant |
| P3 Low Issues | 34 | 🟢 Minor |
| Est. Effort to Fix Code | ~250 hours | — |
| Est. Effort to Fix Docs | ~120 hours | — |

---

## 1. MODULE HEALTH HEATMAP

Accuracy scores based on weighted average: function existence (40%), signature correctness (30%), implementation accuracy (20%), example validity (10%).

| Module | Accuracy | Health | Critical Issues | P1 Issues | Notes |
|--------|----------|--------|----------------|-----------|-------|
| **Ecology** | 96.0% | 🟢 Excellent | 0 | 0 | All formulas verified; minor terminology slip (PCA vs PCoA) |
| **Metagenomics** | 96.2% | 🟢 Excellent | 0 | 0 | One duplicate implementation file |
| **Protein** | 92.9% | 🟢 Good | 1 | 0 | UniProt annotation extraction broken (critical) |
| **Phenotype** | 92.0% | 🟢 Good | 0 | 0 | API exposure issues (deep nesting), missing re-exports |
| **GWAS** | 90.5% | 🟡 Moderate | 2 | 3 | CLI missing, ADMIXTURE misrepresented, QC param issues |
| **Agents** | 89.5% | 🟡 Moderate | 0 | 1 | Truncated file, missing Patterns subsections |
| **Life Events** | 85.0% | 🟡 Moderate | 1 | 2 | Workflow integration breaking mismatch |
| **RNA** | 83.0% | 🟡 Moderate | 0 | 2 | Script path documentation outdated |
| **Simulation** | 76.0% | 🟡 Moderate | 1 | 2 | Fold-change bug, false coalescent claim |
| **ML** | 78.0% | 🔴 Poor | 2 | 4 | Deep learning & LLM APIs missing |
| **SingleCell** | 60.0% | 🔴 Poor | 2 | 5 | DE stubs, trajectory API missing, viz 77% incomplete |
| **Networks** | 45.0% | 🔴 Very Poor | 3 | 6 | 9 functions missing, API confusion throughout |
| **Tasks** | 40.0% | 🔴 Very Poor | 4 | 8 | 40% script paths wrong, 60% flags mismatched |
| **Cross-Code** | N/A | 🔴 Critical | 3,141 violations | — | 1,442 AttributeErrors, 1,263 ImportErrors, 429 SyntaxErrors |

**Heatmap Key:** 🟢≥90% (Healthy) | 🟡70–90% (Moderate) | 🔴<70% (Poor) | 🔴 Critical (Systemic)

---

## 2. TOP 10 CRITICAL ISSUES — RANKED BY IMPACT

### P0 CRITICAL — Must Fix Before Release

| Rank | ID | Module | Title | Severity | Impact |
|------|----|--------|-------|----------|--------|
| 1 | C3 | SingleCell | Differential Expression — Statistical Tests Are Stubs | Critical | Returns random p-values; ALL DE results scientifically invalid |
| 2 | C5 | Simulation | RNA-seq Fold-Change Sign Bug | Critical | Inverts downregulated genes; invalidates simulation benchmarks |
| 3 | C4 | Protein | UniProt Annotation Extraction — Always Returns Empty | Critical | GO terms & keywords always empty; function broken |
| 4 | C6 | GWAS | CLI Command `gwas run` Not Implemented | Critical | Documented command returns error; users cannot follow workflow |
| 7 | C9 | SingleCell | Trajectory Analysis — API Doesn't Exist | Critical | Entire trajectory module missing; docs are fiction |
| 8 | C10 | SingleCell | Visualization — Only 3 of 13 Functions Implemented | Critical | 77% of promised plotting functions unavailable |
| 11 | C1 | ML | Deep Learning — Core API Missing | Critical | embed_dna_sequences, embed_sequences, fine_tune don't exist |
| 12 | C2 | ML | LLM Integration — Domain Methods Missing | Critical | query_biological, interpret_gwas, etc. not implemented |
| 15 | C8 | Simulation | False Coalescent Claim | Critical | Documentation claims sophisticated model; reality is toy simulation |
| 6 | C7 | GWAS | QC Parameter `min_call_rate` Silently Ignored | High | Config mapped but never applied; results include low-call variants |

**Top 3 fastest fixes with highest impact:**
1. Fix fold-change sign bug (C5) — 30 minutes
2. Fix UniProt extraction (C4) — 2 hours
3. Implement statistical tests (C3) — 4 hours (scipy integration)

---

## 3. RISK ASSESSMENT BY DOMAIN

### 🔴 EXTREME RISK — Data Integrity & Scientific Validity
- **SingleCell DE (C3)**: Random p-values produce false discoveries; could invalidate published research
- **Simulation Fold-Change (C5)**: Systematically inverts differential expression direction
- **UniProt Extraction (C4)**: Always returns empty annotations; silent failure mode
- **Risk**: Direct compromise of research reproducibility and scientific integrity

### 🔴 CRITICAL RISK — Functionality & User Trust
- **Deep Learning (C1) & LLM (C2)**: Entire documented features missing
- **GWAS CLI (C6)**: Users cannot execute documented workflows
- **Trajectory Analysis (C9)**: Major documented module completely missing
- **SingleCell Viz (C10)**: 77% of visualizations unavailable
- **Risk**: Users lose confidence; documentation becomes unreliable reference

### 🟠 HIGH RISK — Configuration Drift & Silent Failures
- **GWAS QC params (C7)**: `min_call_rate`, `min_qual`, `exclude_indels` documented but ignored
- **Parameter mismatches**: `max_missing` default 0.05 (docs) vs 0.1 (code), `hwe_pval` vs `min_hwe_p`
- **Risk**: Results silently differ from user expectations; reproducibility compromised

### 🟡 MODERATE RISK — Developer Experience & Maintenance
- **3,141 cross-code violations**: Documentation code examples cannot be executed
- **302 broken links**: Navigation and reference failures across 123 files
- **Tasks module broken examples**: 40% of script paths incorrect
- **Risk**: Onboarding friction; increased support burden; time wasted troubleshooting

### 🟢 LOW RISK — Cosmetic & Polish
- Terminology slips (PCA vs PCoA)
- Missing `__init__` re-exports
- Broader placeholder URLs
- **Risk**: Minor confusion; professional appearance compromised

---

## 4. IMMEDIATE ACTION ITEMS (NEXT 48 HOURS)

### Phase 1 — Scientific Integrity Fixes *(Priority: MUST COMPLETE)*

| # | Task | Module | Est. Time | Owner | Dependencies |
|---|------|--------|-----------|-------|--------------|
| 1.1 | Fix fold-change sign in `simulate_differential_expression()` | Simulation | 0.5h | Eng | None |
| 1.2 | Fix `fetch_uniprot_record()` to extract GO terms & keywords | Protein | 2h | Eng | None |
| 1.3 | Implement `_wilcoxon_rank_sum()` and `_welch_t_test()` using scipy.stats | SingleCell DE | 4h | Eng | scipy installed |
| 1.4 | CRITICAL DECISION: Implement `gwas run` subcommand OR remove docs | GWAS CLI | 2–8h | Eng/PM | Decision needed |

### Phase 2 — Broken Documentation Examples *(Priority: HIGH)*

| # | Task | Module | Est. Time | Owner |
|---|------|--------|-----------|-------|
| 2.1 | Audit and fix `docs/tasks/*.md` broken examples (8 files) | Tasks | 8h | Docs |
| 2.2 | Fix major broken examples: `convert.vcf_to_bed()`, `memoize_disk`, `parallel_map` | Core | 4h | Docs/Eng |
| 2.3 | Correct GWAS visualization examples to match actual API | GWAS | 2h | Docs |
| 2.4 | Fix phenotype import examples to show correct paths | Phenotype | 2h | Docs |

### Phase 3 — Silent Configuration Failures *(Priority: HIGH)*

| # | Task | Module | Est. Time | Owner |
|---|------|--------|-----------|-------|
| 3.1 | Implement `min_call_rate` filter in GWAS QC OR remove from config/docs | GWAS | 2h | Eng |
| 3.2 | Implement `min_qual` and `exclude_indels` filters OR remove from config/docs | GWAS | 3h | Eng |
| 3.3 | Align `max_missing` default: change code to 0.05 OR update docs to 0.1 | GWAS | 0.5h | Eng/Docs |

**48-Hour Total Minimum (Code Fixes Only):** 26 hours (1.1–1.3, 3.1–3.3)
**48-Hour Total Including Task Fixes:** 36 hours

---

## 5. 30-60-90 DAY ROADMAP

### 🎯 PRIORITY 0 — CRITICAL BLOCKERS (Weeks 1–2) — ~89h impl OR ~36h docs

**GOAL:** Unbreak core functionality and restore documentation reliability

| Week | Sprint Focus | Tasks | Effort (Impl) | Effort (Docs) |
|------|--------------|-------|---------------|---------------|
| W1–W2 | SingleCell Scientific Validity | DE statistical tests (C3); Fix 10 missing viz (C10); Implement trajectory API (C9) | 76h | 12h |
| W1 | Simulation Data Integrity | Fold-change fix (C5) | 0.5h | — |
| W1 | Protein API Restoration | UniProt fix (C4) | 2h | — |
| W1–W2 | GWAS CLI & QC | `gwas run` subcommand (C6); `min_call_rate` (C7); `min_qual`/`exclude_indels` | 13h | 4h |
| W2 | Deep Learning & LLM Decision | Review scope; decide implement vs scale back docs (C1, C2) | Planning | 8h |

**Decision Required:**
- **Option A (Implement)**: Deep learning APIs (24h), LLM domain methods (16h) → total 40h
- **Option B (Downscope docs)**: Rewrite deep_learning.md, llm_integration.md → total 8h

---

### 🟡 PRIORITY 1 — HIGH SEVERITY (Sprint 1–2) — ~162h impl OR ~64h docs

**GOAL:** Resolve major API mismatches and implementation gaps

| Week | Sprint Focus | Key Tasks | Effort |
|------|--------------|-----------|--------|
| W3–W4 | Networks Module Overhaul | Implement 9 missing functions OR major docs rewrite (PageRank, spectral, hierarchical, pathway, PPI) | 40h / 40h |
| W3 | AutoML API Alignment | Rewrite automl.md to match model_fn pattern, metric vs scoring | 8h docs |
| W3 | ML Evaluation Docs | Update evaluation.md to match actual function names (`cross_validation_scores`, `bootstrap_validate`) | 6h docs |
| W3–W4 | Tasks Module Cleanup | Audit & rewrite all 8 task guides to match actual scripts | 24h docs |
| W4 | MCP Completion Decision | Complete server (60h) OR mark as experimental (2h) | 60h / 2h |
| W3–W5 | Phenotype Re-exports | Fix all 8 subpackage `__init__.py` to re-export classes; move misplaced `mappings.py` | 9h impl |
| W5 | MLM Review Correction | Fix `review.md` GWAS MLM status (per audit: is implemented) | 1h |
| W3 | Simulation Missing Features | Document `simulate_admixture`, `simulate_selection` in popgen.md | 2h |

**Total Time:** ~162h implementation OR ~64h documentation fix

---

### 🟢 PRIORITY 2 — MEDIUM / POLISH (Sprint 3+) — ~72h additional

**GOAL:** Incremental improvements; lower user friction; address technical debt

| Category | Tasks | Effort |
|----------|-------|--------|
| **Broken Links** | Fix all 288 internal broken links + 14 external | 24h |
| **API Standardization** | Standardize simulation parameter names (maf_min vs min_maf, n_snps vs n_sites) | 4h |
| **Missing Parameters** | Add `gc_content`, `within_pop_diversity`, `recovery_generations`, `growth_rate` or remove from docs | 2–3h each |
| **Protein Enhanced Docs** | Expand proteomes.md to cover 5 undocumented functions; add accession regex for A0A* | 3.5h |
| **Phenotype Metrics** | Document `analysis/statistical.py` functions; add top-level convenience imports | 9h |
| **LD Decay Fix** | Implement proper LD decay in `generate_linkage_disequilibrium_data()` | 8h |
| **Other Simulation** | Forward-time Wright-Fisher / proper coalescent (long-term R&D) | — |

---

### ⚪ PRIORITY 3 — LONG-TERM / ARCHITECTURAL (Q3–Q4)

| # | Initiative | Rationale |
|---|------------|-----------|
| 3.1 | Replace simplified population genetics with forward-time Wright-Fisher or coalescent | Research-grade simulation quality |
| 3.2 | Create auto-generated Sphinx API reference from docstrings | Prevent drift |
| 3.3 | Add comprehensive type hints documentation pages for all modules | Completeness |
| 3.4 | Parallelize batch AlphaFold downloads with ThreadPoolExecutor | Performance |
| 3.5 | Standardize network error handling (raise vs return []) | Consistency |
| 3.6 | Write integration tests validating simulation statistical properties | Quality assurance |

---

## 6. INVESTMENT SCENARIOS

### Scenario A: Full Code Implementation (Maximal Scope) — ~250 engineering hours
**Pros:** Documentation becomes truth; users get all promised features.
**Cons:** Significant engineering investment; may exceed project scope.

| Feature Set | Hours | Modules Improved |
|-------------|-------|------------------|
| Deep Learning embedding + fine-tuning APIs | 24 | ML |
| LLM domain methods (7 functions) | 16 | ML |
| SingleCell DE statistical tests | 4 | SingleCell |
| SingleCell Viz (10 missing plots) | 40 | SingleCell |
| SingleCell Trajectory API | 32 | SingleCell |
| GWAS QC filters (3 missing) | 5 | GWAS |
| GWAS `gwas run` CLI | 8 | GWAS |
| Simulation fixes (fold-change + coalescent/LD/Fst) | 24.5 | Simulation |
| Networks (9 missing functions) | 40 | Networks |
| MCP server completion | 60 | MCP |
| **Total** | **253.5h** | — |

---

### Scenario B: Documentation Correction (Practical Path) — ~120 documentation hours
**Pros:** Faster, lower risk; aligns docs with reality.
**Cons:** Removes features from documentation scope; users lose expectations of these capabilities.

| Task | Hours | Modules Improved |
|------|-------|------------------|
| Rewrite deep_learning.md to match low-level API | 4 | ML |
| Rewrite llm_integration.md to basic wrapper | 4 | ML |
| Rewrite trajectory.md to match actual `dpt_trajectory`/`paga_trajectory` | 6 | SingleCell |
| Rewrite visualization.md to list only 3 implemented plots | 6 | SingleCell |
| Tasks: Audit & rewrite all 8 task guides to match scripts | 24 | Tasks |
| MCP: Mark as experimental draft (40% complete) | 2 | MCP |
| Fix 288 broken internal links | 24 | All |
| Fix broken examples catalog (60+ items) | 30 | All |
| GWAS QC param docs alignment (min_call_rate, min_qual, exclude_indels, max_missing) | 4.5 | GWAS |
| Networks major docs rewrite for 9 missing functions | 40 | Networks |
| Other P1/P2 fixes (AutoML, Phenotype re-exports, Simulation params) | ~30 | Mixed |
| **Total** | **~174.5h** | — |

**Recommendation (Hybrid Approach):**
1. **Fix Code (Scenario A subset, 26h minimum):** C3, C4, C5, C7, C10 fix implementation where it's <2h each
2. **Scale Back Docs (Scenario B subset):** Deep learning, LLM, trajectory, MCP → downscope
3. **Hybrid Total:** ~50h code fixes + ~100h docs fixes = **150h total** over 6–8 weeks

---

## 7. KEY FINDINGS SUMMARY

### What the Numbers Reveal

| Aspect | Finding | Implication |
|--------|---------|-------------|
| **Overall health** | 72.8% composite accuracy | Documentation moderately reliable but significant gaps |
| **Best modules** | Ecology (96%), Metagenomics (96.2%), Protein (92.9%), Phenotype (92%) | These are production-ready; use as reference patterns |
| **Critical failings** | SingleCell (60%), Networks (45%), Tasks (40%), Cross-Code (3,141 violations) | Major rework needed before release |
| **Scientific integrity** | DE stubs returning random p-values (C3), fold-change inversion (C5) | Research published using these tools could be invalid |
| **Documentation-execution gap** | 3,141 cross-code violations across docs; 60+ broken examples | Code examples can't be trusted; users encounter immediate failures |
| **Maintenance debt** | 288 broken internal links; pervasive `AttributeError`s from missing re-exports | Eroding trust; expensive to fix incrementally |

### What Needs to Change

**Before Release (Non-negotiable):**
- Fix DE statistical tests (C3) or clearly mark as stubs
- Fix fold-change sign (C5)
- Fix UniProt extraction (C4) or mark non-functional
- Fix GWAS CLI execution path (C6): either implement `gwas run` or rewrite all workflow docs
- Address silent QC filter skips (C7, C9, C10)

**Short-term (0–3 months):**
- Decide Deep Learning & LLM scope: complete implementation or scale back promises
- Overhaul Tasks documentation (systematic)
- Complete or de-scope MCP server
- Resolve Networks module mismatch (implementation or documentation)

**Long-term (3–6 months):**
- Fix cross-code violations (re-export patterns, module structure)
- Standardize configuration vs API parameter mapping
- Implement population genetics upgrades (coalescent, proper LD/Fst)
- Establish automated documentation testing to prevent future drift

---

## 8. RECOMMENDATION FOR LEADERSHIP

**Strategic Decision Required:** Choose **implementation track** OR **documentation alignment track** — do not attempt both simultaneously unless resources are increased.

### Recommended Path (Balanced — 150h, 6–8 weeks)

**Week 1–2 — Scientific Integrity Sprint (26h eng + 10h docs):**
1. Fix all scientific accuracy bugs (C3, C4, C5) — 6.5h
2. Fix silent config skips (C7, C9, C10 minimal) — 5.5h
3. Fix GWAS CLI or rewrite docs (Decision required, ~2h docs OR 8h impl)
4. Document clear warnings for unimplemented features (Deep Learning, LLM, Trajectory) — 10h

**Week 3–4 — Critical Documentation Repairs (8h eng + 64h docs):**
1. Rewrite Tasks guides (8 files) — 24h docs
2. Rewrite Networks docs for missing API — 40h docs
3. Fix AutoML, Phenotype, Simulation parameter mismatches — 8h docs
4. Fix 60+ broken examples catalog — 30h docs

**Week 5–6 — Link Hygiene & Post-Mortem (24h docs + planning):**
1. Fix all 302 broken internal/external links — 24h docs
2. Post-release retrospective: plan long-term fixes (simulation, MCP, cross-code violations)

**Gate:** Hold release until all P0 critical issues addressed and documentation examples validated executable.

---

## APPENDIX: DATA SOURCES

This dashboard synthesized findings from:
1. **DOCUMENTATION_AUDIT_REPORT.md** — Master cross-module audit (1,154 lines)
2. **17+ module validation reports** — Per-domain accuracy analysis (ecology, metagenomics, protein, phenotype, gwas, agents, life_events, rna, simulation, ml, singlecell, networks, tasks, etc.)
3. **cross_code_verification_report.md** — 3,141 import/attribute/syntax violations
4. **LINK_VALIDATION_REPORT.md** — 302 broken internal/external links across 123 files
5. **Domain-specific issue catalogs** — BROKEN_EXAMPLES.md, SIGNATURE_MISMATCHES.md, etc.

All data current as of report date: April 29, 2026.

---

**END OF EXECUTIVE SUMMARY DASHBOARD**
