# Plan: Documentation Accuracy & Completeness Audit Fix

> **ROUND-5 STATUS PASS (2026-09-01, lane R5-DOCACC)** — every item re-verified against current code:
> 1. DONE (import path already `metainformant.math.population_genetics` in docs/math/selection.md).
> 2. DONE-obsolete (docs/comprehensive_review_report.md, documentation_review_report.md, documentation_audit.md no longer exist at docs/ root — archived under docs/archive/historical_reports/).
> 3. DONE (stale banner present in docs/math/REVIEW_AND_IMPROVEMENT_PLAN.md).
> 4. DONE-obsolete (`setup_uv.sh` has zero matches in docs/, scripts/, CLAUDE.md, README.md).
> 5. VERIFIED-CURRENT (src/metainformant/rna/engine/progress_tracker.py still exists and is referenced by tests/rna/test_rna_amalgkit_current.py — docs/LINUX_TRANSFER.md note stays accurate).
> 6. DONE (`cloud/` present in CLAUDE.md module structure, line 97).
> 7. DONE (src/metainformant/mcp/__init__.py exists; package is proper).
> 8. CONFIRMED-NO-FIX (GWAS docs accurately describe real limitations).
> 9. CONFIRMED-NO-FIX (np.random placeholders are example/demo data only).
> Module-count check: `ls src/metainformant/` = 28 module dirs. Plan COMPLETE — remaining work tracked in lane report.

## Context

Full audit of `docs/` revealed the documentation is overall excellent (415+ files, 30+ subdirectories, consistent structure). However, several issues need fixing: wrong import paths, stale review reports referencing nonexistent files, outdated improvement plans describing already-fixed code, legacy script references in AGENTS.md, a missing module in CLAUDE.md, and a backward-compat mention that should be cleaned up.

---

## Issues & Fixes

### 1. Wrong import path in `docs/math/selection.md`
- **Line 6**: `from metainformant.math import selection` is WRONG
- **Correct**: `from metainformant.math.population_genetics import selection`
- The functions `kin_selection_response` and `multilevel_selection_decomposition` exist in `src/metainformant/math/population_genetics/selection.py`

### 2. Stale review reports reference nonexistent file
- `docs/comprehensive_review_report.md:4` references `docs/comprehensive_docs_review.md` (March 2026) -- this file does NOT exist
- **Fix**: Remove the broken cross-reference, keep the stale warning but remove the "See X for latest" pointer

### 3. `docs/math/REVIEW_AND_IMPROVEMENT_PLAN.md` -- outdated findings
- Claims `fisher_exact_test()` returns placeholder p-value (1.0) -- **this is fixed** in the actual code (uses chi-square approximation now)
- Claims README.md references non-existent `kin_selection.py` / `multilevel_selection.py` -- the src README.md is actually correct now (doesn't reference these)
- **Fix**: Add a stale warning banner similar to the other review reports, noting it was written when these issues existed but many have since been resolved

### 4. `scripts/package/AGENTS.md` -- lists deprecated `setup_uv.sh`
- Line 11: `- setup_uv.sh - Initialize UV environment`
- The script still exists (redirects to setup.sh) but documenting it as a key script is misleading
- **Fix**: Mark it as deprecated redirect, or remove it and note `setup.sh` is the canonical entry point

### 5. `docs/LINUX_TRANSFER.md:65` -- backward compat mention
- `progress_tracker.py` marked as "Deprecated (JSON-based, kept for backward compat)"
- The file still exists. This is accurate documentation of current state.
- **Fix**: Verify if progress_tracker.py is still imported anywhere. If not, note it can be removed. If still used, leave as-is.

### 6. CLAUDE.md missing `cloud/` module
- `src/metainformant/cloud/` exists with `__init__.py` but is not listed in the CLAUDE.md architecture diagram
- **Fix**: Add ` cloud/ # GCP deployment, Docker pipelines` to the module structure

### 7. `src/metainformant/mcp/` missing `__init__.py`
- Directory exists with `tools/` subdirectory but no `__init__.py` -- not a proper Python package
- **Fix**: Create `__init__.py` to make it a proper package, or note in docs that it's a tools directory not a module

### 8. Placeholder mentions in GWAS docs are accurate
- `docs/gwas/index.md:196`, `docs/gwas/ASSESSMENT.md`, `docs/gwas/workflow.md` -- these document real limitations (variant download is a placeholder)
- **No fix needed** -- these are accurate documentation of current state with workarounds noted

### 9. `docs/visualization/ontology.md:241` and `docs/visualization/epigenome.md:185`
- Use `np.random.rand()` as placeholder data in code examples
- **No fix needed** -- these are example/demo data in documentation, not placeholder implementations

---

## Files to Modify

| File | Change |
|------|--------|
| `docs/math/selection.md` | Fix import path |
| `docs/comprehensive_review_report.md` | Remove broken cross-ref to nonexistent file |
| `docs/math/REVIEW_AND_IMPROVEMENT_PLAN.md` | Add stale banner, note fisher_exact_test is fixed |
| `scripts/package/AGENTS.md` | Mark setup_uv.sh as deprecated redirect |
| `CLAUDE.md` (project root) | Add cloud/ to module structure |
| `src/metainformant/mcp/__init__.py` | Create if mcp should be a module (verify first) |

## Files to Potentially Remove (confirm with user)

| File | Reason |
|------|--------|
| `docs/comprehensive_review_report.md` | Stale (Jan 2025), references nonexistent successor |
| `docs/documentation_review_report.md` | Stale (Dec 2025), issues already fixed |
| `docs/documentation_audit.md` | Stale (Dec 2025), fixes already applied |

---

## Verification

1. `grep -r "from metainformant.math import selection" docs/` -- should return 0 matches after fix
2. `grep -r "setup_uv.sh" docs/ scripts/` -- only appears in AGENTS.md with deprecation note
3. `grep -r "comprehensive_docs_review.md" docs/` -- should return 0 matches after fix
4. Confirm CLAUDE.md module list matches `ls src/metainformant/` directories with `__init__.py`
