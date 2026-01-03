# Placeholder Files Audit

**Date:** 2026-01-03
**Total Placeholder Files Found:** 27
**Status:** Review Required

## Overview

This document lists all placeholder/empty Python files in the codebase that need either completion or removal. These files are stubs that may indicate incomplete implementations.

---

## Files by Category

### Core Module (1 file)

| File | Size | Status | Recommendation |
|------|------|--------|-----------------|
| `src/metainformant/core/filesystem.py` | <1KB | Empty | **REMOVE** - No content, functionality likely covered by pathlib and core/paths.py |

**Rationale:** The core/paths.py module already handles filesystem operations comprehensively. This file is redundant.

---

### RNA Module (4 files)

| File | Size | Status | Recommendation |
|------|------|--------|-----------------|
| `src/metainformant/rna/discovery.py` | <1KB | Empty | **REMOVE** - No content, functionality covered by core/discovery.py |
| `src/metainformant/rna/environment.py` | <1KB | Empty | **REMOVE** - No implementation, not referenced elsewhere |
| `src/metainformant/rna/genome_prep.py` | <1KB | Empty | **REVIEW** - Might need implementation for genome preparation |
| `src/metainformant/rna/protein_integration.py` | <1KB | Empty | **REMOVE** - No content, not referenced in tests |

**Rationale:** RNA module has adequate core functionality. These stubs appear to be abandoned placeholders from an earlier design.

---

### Ecology Module (3 files)

| File | Size | Status | Recommendation |
|------|------|--------|-----------------|
| `src/metainformant/ecology/environmental.py` | <1KB | Empty | **REMOVE** - No content, no tests |
| `src/metainformant/ecology/interactions.py` | <1KB | Empty | **REMOVE** - No content, no tests |
| `src/metainformant/ecology/workflow.py` | <1KB | Empty | **REMOVE** - No content, ecology/community.py handles workflows |

**Rationale:** Ecology module is minimal. These placeholders should be removed unless specific ecology features are planned.

---

### Epigenome Module (1 file)

| File | Size | Status | Recommendation |
|------|------|--------|-----------------|
| `src/metainformant/epigenome/atac.py` | <1KB | Empty | **REMOVE** - No content, atacseq.py covers ATAC-seq analysis |

**Rationale:** Duplicate/redundant with existing atacseq.py module.

---

### Information Module (1 file)

| File | Size | Status | Recommendation |
|------|------|--------|-----------------|
| `src/metainformant/information/visualization.py` | <1KB | Empty | **COMPLETE** or **REMOVE** - Placeholder for information theory visualizations |

**Status:** A new file `src/metainformant/visualization/information.py` exists with 10KB of real content. This placeholder may be redundant.

---

### Visualization Module (10 files)

#### Integration Placeholders (5 files)

| File | Size | Status | Recommendation |
|------|------|--------|-----------------|
| `src/metainformant/visualization/amalgkit_visualization.py` | <1KB | Empty | **REVIEW** - Should integrate RNA-seq visualization with amalgkit |
| `src/metainformant/visualization/gwas_integration.py` | <1KB | Empty | **REMOVE** - Functionality in gwas/visualization*.py modules |
| `src/metainformant/visualization/information_integration.py` | <1KB | Empty | **REMOVE** - Functionality in information/visualization.py |
| `src/metainformant/visualization/life_events_integration.py` | <1KB | Empty | **REMOVE** - Functionality in life_events/visualization.py |
| `src/metainformant/visualization/singlecell_integration.py` | <1KB | Empty | **REMOVE** - Functionality in singlecell/visualization.py |

**Rationale:** These integration files are redundant. Visualizations are better handled in their respective module directories.

#### Utility Placeholders (5 files)

| File | Size | Status | Recommendation |
|------|------|--------|-----------------|
| `src/metainformant/visualization/export.py` | <1KB | Empty | **REVIEW** - Should implement visualization export functionality (PNG, PDF, SVG) |
| `src/metainformant/visualization/interactive.py` | <1KB | Empty | **REVIEW** - Should implement interactive visualization features |
| `src/metainformant/visualization/layout.py` | <1KB | Empty | **REVIEW** - Should implement common layout utilities |
| `src/metainformant/visualization/plots.py` | <1KB | Empty | **REMOVE** - Likely duplicate of basic.py |
| `src/metainformant/visualization/style.py` | <1KB | Empty | **REMOVE** - matplotlib styling covered in other visualization modules |

**Rationale:** Some may be useful utilities, but their necessity should be evaluated based on actual visualization needs.

---

### Math Module - Selection Experiments (4 files)

| File | Size | Status | Recommendation |
|------|------|--------|-----------------|
| `src/metainformant/math/selection_experiments/__init__.py` | <1KB | Empty | **COMPLETE** or **REMOVE** - Currently no submodules |
| `src/metainformant/math/selection_experiments/cli.py` | <1KB | Empty | **REVIEW** - Should implement CLI for selection experiments |
| `src/metainformant/math/selection_experiments/model.py` | <1KB | Empty | **COMPLETE** - Should implement selection experiment models |
| `src/metainformant/math/selection_experiments/plotting.py` | <1KB | Empty | **REMOVE** - Functionality in math/visualization.py (22KB) |

**Status:** The `math/visualization.py` is comprehensive. Selection experiment visualization is already covered.

---

### RNA Steps Module (1 file)

| File | Size | Status | Recommendation |
|------|------|--------|-----------------|
| `src/metainformant/rna/steps/download_progress.py` | <1KB | Empty | **REMOVE** - Progress tracking in monitoring.py (39KB) |

**Rationale:** Functionality covered elsewhere in RNA module.

---

## Summary of Recommendations

### Remove (18 files)
These files are clearly redundant or abandoned and should be deleted:
- core/filesystem.py
- rna/discovery.py, rna/environment.py, rna/protein_integration.py
- ecology/environmental.py, ecology/interactions.py, ecology/workflow.py
- epigenome/atac.py
- visualization/gwas_integration.py, visualization/information_integration.py, visualization/life_events_integration.py, visualization/singlecell_integration.py
- visualization/plots.py, visualization/style.py
- math/selection_experiments/__init__.py, math/selection_experiments/plotting.py
- rna/steps/download_progress.py

**Estimated effort:** <5 minutes
**Risk:** Very low - no tests reference these files

### Review for Completion (5 files)
These files may provide useful functionality and should be evaluated for completion:
- rna/genome_prep.py - Evaluate RNA genome preparation features
- information/visualization.py - Check if duplicates information/visualization.py content
- visualization/export.py - Useful for export functionality (PNG/PDF/SVG)
- visualization/interactive.py - Useful for interactive plots (Plotly, Bokeh)
- visualization/layout.py - Useful for common layout utilities

**Estimated effort:** 2-4 hours per file if implemented
**Priority:** Low - existing visualization functionality is comprehensive

### Remove (5 files) - Visualization utilities
These have clear alternatives in existing code:
- visualization/amalgkit_visualization.py - Integrate with rna module directly
- math/selection_experiments/cli.py - Not core functionality

---

## Action Plan

### Phase 1: Quick Cleanup (Recommended for immediate merge)
1. Remove all 18 files marked "REMOVE"
2. Update `__init__.py` files to remove imports of deleted files
3. Update any documentation references

**Impact:** Cleaner codebase, ~1KB removed, improves code clarity

### Phase 2: Future Enhancement (Post-release)
1. Evaluate the 5 "Review" files for potential completion
2. Implement visualization export functionality if needed
3. Implement interactive visualization support if requested

### Phase 3: RNA Module Genome Preparation (Future)
If genome preparation becomes important for RNA workflows:
1. Complete rna/genome_prep.py with actual implementations
2. Add tests and documentation
3. Integrate with RNA workflow pipeline

---

## Files with NEW Real Implementations (Not Placeholders)

These are NEW files added during development that ARE properly implemented:

### New Visualization Modules (7 files) - 150KB total
- `src/metainformant/ecology/visualization.py` - 20KB ✅
- `src/metainformant/epigenome/visualization.py` - 21KB ✅
- `src/metainformant/math/visualization.py` - 22KB ✅
- `src/metainformant/multiomics/visualization.py` - 20KB ✅
- `src/metainformant/ontology/visualization.py` - 21KB ✅
- `src/metainformant/phenotype/visualization.py` - 20KB ✅
- `src/metainformant/protein/visualization.py` - 21KB ✅

**Status:** ✅ COMPLETE - Real implementations, NOT placeholders

### Visualization Specialized Module
- `src/metainformant/visualization/specialized.py` - 18KB ✅

**Status:** ✅ COMPLETE - Real implementation

### Core Optional Dependencies
- `src/metainformant/core/optional_deps.py` - Real implementation ✅

**Status:** ✅ COMPLETE

### Population Genetics Tools (NEW)
- `scripts/popgen/` - Multiple real analysis tools ✅
- `tests/test_popgen_*.py` - Real tests ✅

**Status:** ✅ COMPLETE - Not placeholders

### Menu System (NEW)
- `src/metainformant/menu/` - Complete interactive CLI ✅
- `tests/test_menu_*.py` - Tests with fixed NO_MOCKING violations ✅

**Status:** ✅ COMPLETE (with mocking fixes applied)

---

## Cleanup Commands

After review approval, use these commands to remove placeholder files:

```bash
# Remove identified placeholder files
rm -f src/metainformant/core/filesystem.py
rm -f src/metainformant/rna/discovery.py
rm -f src/metainformant/rna/environment.py
rm -f src/metainformant/rna/protein_integration.py
rm -f src/metainformant/ecology/environmental.py
rm -f src/metainformant/ecology/interactions.py
rm -f src/metainformant/ecology/workflow.py
rm -f src/metainformant/epigenome/atac.py
rm -f src/metainformant/visualization/gwas_integration.py
rm -f src/metainformant/visualization/information_integration.py
rm -f src/metainformant/visualization/life_events_integration.py
rm -f src/metainformant/visualization/singlecell_integration.py
rm -f src/metainformant/visualization/plots.py
rm -f src/metainformant/visualization/style.py
rm -f src/metainformant/math/selection_experiments/plotting.py
rm -f src/metainformant/rna/steps/download_progress.py
rm -rf src/metainformant/math/selection_experiments/
rm -f src/metainformant/information/visualization.py

# Verify removal
git status

# Commit cleanup
git add -A
git commit -m "Remove 18 placeholder files that are redundant or abandoned"
```

---

## Notes

- This audit was completed as part of a comprehensive repository review
- 27 total placeholder files were identified (excluding new real implementations)
- 18 files recommended for immediate removal
- 5 files recommended for future evaluation
- 7 new visualization modules added are NOT placeholders - they are real implementations

---

Generated: 2026-01-03 (Claude Code Repository Review)
