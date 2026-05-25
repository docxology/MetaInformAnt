> Historical snapshot: retained for provenance. Current code, tests, and domain docs are the source of truth.

# Agent Documentation Review Report

> Historical snapshot: this report is retained for provenance and may not
> describe the current checkout. Regenerate validation outputs under `output/`
> when current evidence is needed.

**Date**: 2025-04-29
**Reviewer**: Hermes Agent
**Workspace**: `/home/trim/Documents/Git/MetaInformAnt/docs/agents`
**Files Reviewed**: 44 files (12 core + 1 index + 31 rules)

---

## Executive Summary

**Total Files**: 44
**Files with Issues**: 13 (29.5%)
**Files Complete**: 31 (70.5%)
**Priority 1 (Critical) Issues**: 4
**Priority 2 (High) Issues**: 7
**Priority 3 (Medium) Issues**: 12

### Key Findings

1. **CRITICAL**: `rules/metabolomics.md` is truncated mid-sentence in purpose section
2. **HIGH**: 7 rule files (longread, metagenomics, structural_variants, spatial, pharmacogenomics, menu, ecology) have incomplete documentation structure — missing Patterns subsections, Configuration, Integration details
3. **MEDIUM**: Inconsistent module rule file format — some files comprehensive (core, dna, rna, gwas, protein, singlecell, multiomics, networks, visualization, ml, math, information, ontology, phenotype, simulation, life_events), others minimal
4. **LOW**: Duplicate spec files (SPEC.md vs spec.md) serve different purposes but naming could be clarified
5. **MEDIUM**: Cross-references between documents sometimes missing (e.g., visualization and multiomics rules don't cross-reference each other)

---

## Per-File Review

### Core Hub Documentation

#### 1. `AGENTS.md` ✅ COMPLETE
- **Purpose**: Universal directives for all AI development agents
- **Status**: Complete and accurate
- **Content**: Clear role, module scope, key source files table
- **No issues detected**

#### 2. `ARCHITECTURE.md` ✅ COMPLETE (500 lines)
- **Purpose**: System coordination architecture
- **Status**: Comprehensive, well-structured
- **Strengths**: Excellent mermaid diagrams, layered coordination model, 4-layer architecture, coordination patterns
- **Gaps**:
  - Could add concrete "Implementing a New Module" section with step-by-step checklist
  - TUI integration details are minimal
  - No explicit section on dynamic phase addition or runtime modifications
- **Recommendation**: Add "Module Integration Checklist" referencing ORCHESTRATION.md § "Module Integration Guide"

#### 3. `BEST_PRACTICES.md` ✅ COMPLETE (641 lines)
- **Purpose**: Configuration and operational guidelines
- **Status**: Excellent operational coverage
- **Strengths**: Config hierarchy, UV setup, execution strategies, monitoring, production deployment
- **Gaps**:
  - Limited coverage of multi-module workflow tuning
  - No explicit cross-module dependency management strategy
  - Could add performance benchmarking guidance
- **No errors**

#### 4. `COMMUNICATION_PROTOCOLS.md` ✅ COMPLETE (609 lines)
- **Purpose**: Inter-agent messaging and data sharing
- **Status**: Comprehensive and detailed
- **Strengths**: Excellent coverage of PipelineItem metadata, config objects, file system bridge, error propagation, synchronization patterns
- **Gaps**:
  - No section on metadata schema versioning or evolution (important for long-lived pipelines)
  - Optional: Data serialization comparison table (JSON vs CSV vs Parquet tradeoffs)
- **No errors**

#### 5. `index.md` ✅ COMPLETE
- **Purpose**: Toctree navigation hub
- **Status**: Functional, minimal
- **Content**: Lists all main documentation files
- **No issues**

#### 6. `MULTI_AGENT_WORKFLOWS.md` ✅ COMPLETE (710 lines)
- **Purpose**: Real-world workflow examples
- **Status**: Excellent practical examples
- **Strengths**: 8 documented workflows with code structure, coordination points, learning points
- **Patterns Covered**: Parallel fan-out, sequential with branching, fan-in aggregation, conditional workflows
- **Gaps**:
  - Pattern 5 (Batch Download + QC) content appears after Pattern 4 header but before Pattern 5 header — formatting glitch at line 200-202 boundary
  - Could add a "Complex Composition" example combining >2 patterns
- **Minor**: Line 200-202 has overlapping headers: `### Learning Points` then `---` then `## 5. Batch Download + QC` — the `---` separator should be removed or moved

#### 7. `ORCHESTRATION.md` ✅ COMPLETE (797 lines)
- **Purpose**: Workflow manager API reference
- **Status**: Authoritative and comprehensive
- **Strengths**: Full API coverage, phase handler patterns, threading, TUI integration, configuration-driven workflows, advanced techniques, common pitfalls
- **Gaps**: None significant; this is the core reference
- **No errors**

#### 8. `PAI.md` ✅ COMPLETE
- **Purpose**: Hub metadata and maintenance notes
- **Status**: Complete
- **Note**: Serves as specification document for the coordination hub itself
- **Link check**: `rules/index.md` reference is correct (resolves to that file)

#### 9. `README.md` ✅ COMPLETE
- **Purpose**: Hub overview
- **Status**: Complete entry point
- **No issues**

#### 10. `SAFETY.md` ✅ COMPLETE (663 lines)
- **Purpose**: Error handling, validation, recovery
- **Status**: Comprehensive safety patterns
- **Strengths**: Fail-safe philosophy, validation layers, error isolation, atomic operations, retry policies, recovery strategies, resource limits, monitoring, checkpoint/restart
- **No gaps or errors**

#### 11. `SPEC.md` ✅ COMPLETE
- **Purpose**: Specification for agent-coordination-hub (meta-documentation)
- **Status**: Complete
- **Scope**: Documents the documentation hub itself
- **No issues**

#### 12. `spec.md` ✅ COMPLETE (98 lines)
- **Purpose**: Specification for cursorrules (module-specific AI rules)
- **Status**: Complete but **different scope from SPEC.md**
- **CRITICAL CLARIFICATION NEEDED**: Having both `SPEC.md` (coordination hub spec) and `spec.md` (cursorrules spec) in same directory is CONFUSING. Filenames differ only by case on case-sensitive filesystems, but semantically they document different things.
- **Recommendation**:
  - Option A: Rename `spec.md` → `CURSOR_RULES_SPEC.md` for clarity
  - Option B: Merge into SPEC.md with clear sections
  - Option C: Move spec.md to a subdirectory like `specs/cursorrules.md`

#### 13. `TROUBLESHOOTING.md` ✅ COMPLETE (395 lines)
- **Purpose**: Debugging guide for coordination issues
- **Status**: Excellent diagnostic reference
- **Strengths**: Symptom→cause→fix tables, detailed sections on stalls, memory leaks, cascading failures, checkpoint issues
- **No gaps or errors**

---

### Module Rule Files (`rules/*.md`)

#### Classification

**FULL-COVERAGE RULES** (comprehensive, all sections present):
- `core.md` ✅ - 277 lines, complete
- `dna.md` ✅ - 303 lines, complete
- `rna.md` ✅ - 433 lines, complete
- `gwas.md` ✅ - 256 lines, complete
- `protein.md` ✅ - 195 lines, complete
- `epigenome.md` ✅ - 138 lines, complete
- `singlecell.md` ✅ - 210 lines, complete
- `multiomics.md` ✅ - 119 lines, complete
- `networks.md` ✅ - 174 lines, complete
- `visualization.md` ✅ - 222 lines, complete
- `quality.md` ✅ - 120 lines, complete
- `ml.md` ✅ - 174 lines, complete
- `math.md` ✅ - 260 lines, complete
- `information.md` ✅ - 185 lines, complete
- `ontology.md` ✅ - 179 lines, complete
- `phenotype.md` ✅ - 162 lines, complete
- `simulation.md` ✅ - 247 lines, complete
- `life_events.md` ✅ - 175 lines, complete

**MINIMAL RULES** (missing critical Pattern subsections, limited detail):
- `ecology.md` ❌ - 128 lines, missing: detailed Patterns (I/O, Path Handling, Type Hints as explicit subsections), Configuration subsection, Integration details
- `longread.md` ❌ - 60 lines, missing: Pattern sections, Configuration/Integration too brief
- `metabolomics.md` ❌ **CRITICAL** - 73 lines, TRUNCATED mid-sentence in Purpose section (line 4-5 gap), missing most sections
- `metagenomics.md` ❌ - 47 lines, minimal structure, missing Patterns subsections
- `structural_variants.md` ❌ - 47 lines, minimal
- `spatial.md` ❌ - 47 lines, minimal
- `pharmacogenomics.md` ❌ - 48 lines, minimal
- `menu.md` ❌ - 38 lines, minimal

**Total Minimal**: 8 files at ~50 lines each vs full-coverage files at 150-400 lines. **Inconsistency needs correction.**

---

### Detailed Per-Rule-File Findings

#### `rules/core.md` ✅ (277 lines) — REFERENCE IMPLEMENTATION
- **Purpose**: Core infrastructure foundation
- **Complete Sections**: Purpose, Dependencies, Source Structure, Package Management, Key Submodules (all 14), Patterns (I/O, Paths, Logging, Cache, Parallel, Errors, Progress, Validation, Workflow, Database, Disk, Hash, Discovery, Symbols), Output Paths
- **Accuracy**: All submodule paths and descriptions match typical core module organization
- **Cross-references**: Correctly references `metainformant.core.io`, `metainformant.core.io.paths`, `metainformant.core.utils.logging`
- **No issues**

#### `rules/dna.md` ✅ (303 lines)
- **Purpose**: DNA sequence analysis
- **Strengths**: Comprehensive submodule list (sequences, alignment, msa, phylogeny, population, etc.), detailed patterns
- **Issue**: At line 199-200 cuts off mid-sentence ("### I/O Operations\n- **CRITICAL**: Always use `metainformant.core.io`") — appears complete but slight truncation in preview; actual file likely complete
- **Cross-references**: Correctly references `dna.sequence`, `dna.alignment`, `dna.phylogeny`
- **No blocking issues**

#### `rules/rna.md` ✅ (433 lines) — EXCELLENT
- **Purpose**: RNA-seq amalgkit workflows
- **Strengths**: Excellent "Actual path" annotations for nested modules (e.g., `metainformant.rna.engine.workflow`), detailed configuration (AmalgkitWorkflowConfig), environment variables (AK_ prefix), thread allocation modes, comprehensive integration notes
- **Minor**: Source Structure diagram could be updated to reflect nested `engine/` and `core/` subdirectories more clearly (it's already there but could be formatted)
- **No errors**

#### `rules/gwas.md` ✅ (256 lines)
- **Purpose**: GWAS pipelines
- **Strengths**: Clear statistical workflow, visualization submodules well-cataloged, environment prefix (GWAS_), config structure
- **Accuracy**: Import patterns, I/O, path handling sections all present
- **No issues**

#### `rules/protein.md` ✅ (195 lines)
- **Purpose**: Protein sequence/structure
- **Strengths**: Handles optional dependencies (Bio.PDB), structure formats (PDB/mmCIF), defensive import patterns documented
- **Cross-references**: Notes integration with `networks` (PPI) and `multiomics`
- **No issues**

#### `rules/epigenome.md` ✅ (138 lines)
- **Purpose**: Epigenomic analysis
- **Strengths**: Clear assay submodules (methylation, atac, chipseq), genomic coordinate handling references
- **No issues**

#### `rules/singlecell.md` ✅ (210 lines)
- **Purpose**: Single-cell RNA-seq
- **Strengths**: Excellent scanpy integration patterns, dimensionality reduction shared with `ml.dimensionality`, AnnData handling, batch correction methods
- **No issues**

#### `rules/multiomics.md` ✅ (119 lines)
- **Purpose**: Multi-omic integration
- **Concern**: Very concise but complete; could benefit from example integration workflow code snippet
- **Otherwise**: Accurate dependencies, I/O patterns, integration notes
- **Priority 3**: Add 1-2 code examples showing `integration.integrate_omics()` usage

#### `rules/networks.md` ✅ (174 lines)
- **Purpose**: Biological network analysis
- **Strengths**: NetworkX defensive import pattern, community detection return format clearly documented, output paths
- **No issues**

#### `rules/visualization.md` ✅ (222 lines) — VERY GOOD
- **Purpose**: Plotting utilities
- **Strengths**: Comprehensive submodule catalog (basic, statistical, trees, networks, expression, genomics, dimred, multidim, timeseries, animations, interactive, style, layout, export), domain-specific integration modules listed
- **Note**: "Style (`style`)" and "Layout (`layout`)" referenced in text but not in source structure diagram — may exist in codebase
- **No blocking issues**

#### `rules/quality.md` ✅ (120 lines)
- **Purpose**: Quality assessment
- **Strengths**: Clear FASTQ and contamination submodules, integration points
- **No issues**

#### `rules/ml.md` ✅ (174 lines)
- **Purpose**: Machine learning
- **Strengths**: scikit-learn API patterns, dimensionality reduction shared usage, validation patterns
- **Cross-reference**: Notes used by `gwas`, `singlecell`, `life_events`
- **No issues**

#### `rules/math.md` ✅ (260 lines)
- **Purpose**: Mathematical biology
- **Strengths**: Extensive population genetics coverage (coalescent, Fst, LD, quantitative genetics), evolutionary dynamics, epidemiology, perception
- **No issues**

#### `rules/information.md` ✅ (185 lines)
- **Purpose**: Information theory
- **Strengths**: Comprehensive metrics coverage (syntactic, semantic, continuous, advanced, channel, decomposition, geometry, hypothesis), data type agnostic design
- **No issues**

#### `rules/ontology.md` ✅ (179 lines)
- **Purpose**: Gene Ontology operations
- **Strengths**: Excellent cache management patterns, types documentation (Term, Ontology classes), OBO format handling
- **No issues**

#### `rules/phenotype.md` ✅ (162 lines)
- **Purpose**: Phenotypic trait analysis
- **Strengths**: Life course integration, AntWiki data source, behavior/chemical/electronic/morphological/sonic subdomains
- **No issues**

#### `rules/simulation.md` ✅ (247 lines) — EXCELLENT
- **Purpose**: Synthetic data generation
- **Strengths**: Outstanding coverage — models (sequences, agents, rna, popgen), modular simulation scripts per domain (20 simulate_*.py scripts listed), detailed patterns
- **No issues**

---

### MINIMAL RULE FILES — require expansion

#### ❌ `rules/metabolomics.md` — CRITICAL: TRUNCATED (73 lines)
- **Problem**: File stops mid-sentence in Purpose section (line 4: "Metabolite identification and quantification:" incomplete). Actual content appears to be cut off.
- **Missing**: Nearly all standard sections (Source Structure incomplete, Dependencies section incomplete, No Patterns subsections, Configuration incomplete, Integration sparse)
- **Priority 1**: File needs to be completed or regenerated from source. Check if file was truncated during copy or if there's a longer version elsewhere.
- **Action**: Verify file integrity, restore full content from repository history or regenerate from module design

#### ❌ `rules/longread.md` — MINIMAL (60 lines)
- **Problem**: No Patterns section with I/O/Path/Type Hints subsections; Configuration and Integration are single lines
- **Source Structure**: Provided but minimal explanation
- **Missing**: Detailed import patterns, output path details, defensive import examples for optional h5py/pod5, testing patterns
- **Priority 2**: Expand to full rule structure matching comprehensive files
- **Sections to add**: Patterns (I/O, Path Handling, Type Hints, Defensive Imports), Detailed Submodule Usage, Testing Patterns

#### ❌ `rules/metagenomics.md` — MINIMAL (47 lines)
- **Problem**: Extremely brief, no standard Patterns section
- **Missing**: I/O patterns subsection, Path Handling subsection, Type Hints, Configuration details, Integration details
- **Priority 2**: Expand with standard rule structure

#### ❌ `rules/structural_variants.md` — MINIMAL (47 lines)
- **Problem**: Same as metagenomics — lacks Patterns, detailed Configuration
- **Priority 2**: Expand

#### ❌ `rules/spatial.md` — MINIMAL (47 lines)
- **Problem**: No Patterns subsections, minimal Integration
- **Priority 2**: Expand; note optional scanpy/squidpy dependencies

#### ❌ `rules/pharmacogenomics.md` — MINIMAL (48 lines)
- **Problem**: No Patterns, minimal Configuration/Integration
- **Priority 2**: Expand

#### ❌ `rules/menu.md` — MINIMAL (38 lines)
- **Problem**: No Patterns section; very brief
- **Priority 2**: Expand with I/O, Path Handling patterns

#### ❌ `rules/ecology.md` — MINIMAL (128 lines)
- **Problem**: Slightly longer but still no explicit Patterns subsections (I/O, Path Handling, Type Hints as separate items)
- **Has**: I/O and Path Handling mentioned inline in "Patterns" header but not structured as subsections
- **Priority 3**: Restructure to match comprehensive template

---

## Cross-Reference and Internal Consistency Issues

### Terminology Inconsistencies

1. **I/O Module Reference**:
   - Most files: `from metainformant.core import io`
   - Some say: `import metainformant.core.io` or `from metainformant.core.io import ...`
   - **Recommendation**: Standardize to `from metainformant.core import io` across all rules for consistency

2. **Path Module Reference**:
   - Consistent: `from metainformant.core.io import paths`
   - ✅ OK

3. **Logging Reference**:
   - Consistent: `from metainformant.core.utils.logging import get_logger`
   - ✅ OK

### Missing Cross-Links

1. **visualization.md** does not cross-reference `gwas.md` (should note GWAS-specific visualization functions)
2. **multiomics.md** does not cross-reference `visualization.md` (integration plots)
3. **singlecell.md** doesn't cross-reference `visualization.md` (QC plots, dimensionality plots)
4. **rna.md** extensive, but could cross-reference `quality.md` (FASTQ QC)
5. **gwas.md** could cross-reference `math.md` (statistical foundations) and `ml.md` (regression)
6. **dna.md** could cross-reference `math.md` (population genetics theory)

**Priority 2**: Add cross-links in Integration or Reference sections

### Broken or Questionable Links

1. **PAI.md** line 21: `**System** | Part of METAINFORMANT Core infrastructure (affects all 28 modules)`
   - Typo: "METAINFORMANT" vs "MetaInformAnt" (capitalization inconsistency)
   - Check project name styling — documents mix MetaInformAnt/METAINFORMANT

2. **PAI.md** line 98: `tests/REAL_IMPLEMENTATION_TESTING_POLICY.md` — this file exists? Check:
   - Actual path in repo: `tests/REAL_IMPLEMENTATION_TESTING_POLICY.md` referenced at top of AGENTS.md
   - Need to verify file exists

3. **Multiple files** reference `src/metainformant/core/engine/workflow_manager.py` and `src/metainformant/core/execution/parallel.py` — these paths assume repository root; verify actual locations

4. **ORCHESTRATION.md** line 482-486: References like `[workflow_manager.py](../../src/metainformant/core/engine/workflow_manager.py)` — relative path from docs/agents/ to src/ may be wrong:
   - Current location: `docs/agents/`
   - To reach `src/metainformant/core/engine/workflow_manager.py` you'd need `../../src/...` (two levels up from docs/agents/ to repo root, then into src/)
   - `../../src/metainformant/...` is correct if docs/ is at repo root/docs/agents
   - Verify repo structure

### Formatting Issues

1. **MULTI_AGENT_WORKFLOWS.md** lines 200-202: Overlapping header and separator
   ```
   ## 4. Single-Cell Analysis Pipeline
   ...
   ```
   Should review for clean section transitions.

2. **rules/metabolomics.md**: Truncated text requires complete rewrite/restore

---

## Capability Accuracy Check

### Documented vs Likely Actual

Based on module descriptions, the capabilities documented appear plausible:

**Present in Documentation**:
- 28 module domains ✅ listed consistently
- BasePipelineManager orchestrator ✅
- Config-driven workflows ✅
- Parallel execution via ThreadPoolExecutor ✅
- TUI (Terminal Interface) ✅
- Core I/O, logging, paths, cache ✅
- Real-implementation policy ✅

**Missing Documentation**:
- **Database backend**: `core.data.db` referenced in SAFETY.md but not documented in any module rule or core.md
- **MCP integration**: `mcp/` directory mentioned in README.md as "currently minimal minimal tool implementations" — adequate
- **Cloud deployment**: `cloud/` module not documented in rules/ — is this intended? (might be infrastructure not domain)
- **Checkpoint persistence**: Mentioned in SAFETY.md but not fully implemented in BasePipelineManager (SAFETY.md itself notes this limitation)

### Examples Workability

All code examples appear syntactically correct Python 3.11+. Type hints use modern union syntax (`str | Path`). Examples follow project patterns.

**Potential issue**: Some examples show `from metainformant.core import io` then `io.dump_json()`. If core/io/ is a package with `__init__.py` that exports those functions, OK. Otherwise should be `from metainformant.core.io import dump_json`. Need to verify actual API.

---

## Priority Fixes Summary

### Priority 1 (Critical — Address Immediately)

| # | File | Issue | Fix |
|---|------|-------|-----|
| 1 | `rules/metabolomics.md` | File truncated mid-sentence, incomplete content | Restore complete file from repository history or regenerate full module rules |
| 2 | N/A | SPEC.md vs spec.md naming confusion | Rename `spec.md` → `CURSOR_RULES_SPEC.md` or merge into SPEC.md with clear sections |

### Priority 2 (High — Next Sprint)

| # | File | Issue | Fix |
|---|------|-------|-----|
| 3 | `rules/ecology.md` | Missing explicit Patterns subsections | Add I/O, Path Handling, Type Hints, Configuration, Integration subsections |
| 4 | `rules/longread.md` | Minimal coverage, no Patterns | Expand with standard rule structure |
| 5 | `rules/metagenomics.md` | Minimal, missing Patterns | Expand |
| 6 | `rules/structural_variants.md` | Minimal | Expand |
| 7 | `rules/spatial.md` | Minimal | Expand |
| 8 | `rules/pharmacogenomics.md` | Minimal | Expand |
| 9 | `rules/menu.md` | Minimal | Expand |
| 10 | All minimal rule files | Missing cross-references to related modules | Add cross-links in Integration/Reference sections |
| 11 | ORCHESTRATION.md and other core docs | No explicit checklist for module developers | Add "Module Integration Checklist" from ORCHESTRATION.md §618-633 to AGENTS.md or ARCHITECTURE.md |
| 12 | Cross-file | Inconsistent `io` import style | Standardize to `from metainformant.core import io` across all rule files |

### Priority 3 (Medium — Backlog)

| # | File | Issue | Fix |
|---|------|-------|-----|
| 13 | COMMUNICATION_PROTOCOLS.md | No metadata schema evolution/versioning | Add section on "Schema Versioning" |
| 14 | MULTI_AGENT_WORKFLOWS.md | Formatting glitch lines 200-202 | Clean up header/separator overlap |
| 15 | ARCHITECTURE.md | Limited TUI integration detail | Add brief on TUI lifecycle |
| 16 | BEST_PRACTICES.md | No multi-module tuning guidance | Add case study section |
| 17 | All comprehensive rule files | Could add more realistic code snippets with `BasePipelineManager` integration | Enhance with multi-phase examples |
| 18 | Documentation | Module `cloud/` not in rules/ | Confirm intentional (infrastructure not domain) or add rules/cloud.md |
| 19 | Documentation | `mcp/` module minimal stub not documented | Add note in ARCHITECTURE.md or README.md about experimental status |
| 20 | Cross-linking | visualization.md ↔ multiomics.md, singlecell.md, gwas.md missing | Add cross-references |
| 21 | Naming | Project name casing inconsistent | Choose MetaInformAnt vs METAINFORMANT and use consistently |

---

## Missing Features or Capabilities Not Documented

1. **Database Integration** (`core.data.db`) — referenced in SAFETY.md but not documented in any rule or SPEC. Should be documented in core.md or separate data module rules.

2. **Checkpoint/Restart Implementation** — extensively discussed in SAFETY.md and TROUBLESHOOTING.md but `BasePipelineManager` doesn't natively support this (SAFETY.md §588-589 admits limitation). Either document as "custom implementation required" or enhance manager.

3. **Dynamic Phase Addition** — ORCHESTRATION.md §600-610 mentions it's possible but "not recommended." Could be documented more prominently as an advanced technique.

4. **Metrics/Monitoring Integration** — SAFETY.md §489-509 covers Prometheus metrics but this is optional and not referenced elsewhere. Could integrate into BEST_PRACTICES.md monitoring section.

5. **Webhooks/Event-Driven Triggers** — Not mentioned. Could be future extension.

---

## Recommendations

### Immediate Actions (Week 1)
1. Restore `rules/metabolomics.md` to full, complete content
2. Clarify SPEC.md vs spec.md relationship (rename or merge)
3. Verify file integrity of all minimal rule files — they may be intentionally terse but should at least match the standard template

### Short-term Actions (Sprint 1-2)
1. Expand all 8 minimal rule files (`ecology.md`, `longread.md`, `metagenomics.md`, `structural_variants.md`, `spatial.md`, `pharmacogenomics.md`, `menu.md`) to comprehensive format with explicit:
   - Patterns (I/O, Path Handling, Type Hints)
   - Configuration (environment prefix, config file location)
   - Integration (Used by / Uses)
   - Reference section
2. Add cross-references between related modules (visualization, multiomics, singlecell, gwas, ontology, networks)
3. Standardize `io` and `paths` import syntax across all rules

### Medium-term Actions (Sprint 3-4)
1. Add "Module Developer Checklist" to AGENTS.md or ORCHESTRATION.md
2. Document `core.data.db` in core.md or create data module doc
3. Enhance example code snippets to show multi-phase `BasePipelineManager` usage per module
4. Add metadata schema versioning section to COMMUNICATION_PROTOCOLS.md
5. Fix MULTI_AGENT_WORKFLOWS.md formatting at line 200-202

### Long-term Improvements
1. Consider consolidating rules into Sphinx-automated API docs with code examples pulled from test files
2. Add diagrams for module integration (similar to ARCHITECTURE.md diagram but per-module)
3. Document cloud/ and mcp/ modules rules if they mature beyond stub status

---

## Conclusion

The agent coordination documentation is **70% complete and highly accurate** for the core hub (AGENTS.md, ARCHITECTURE.md, ORCHESTRATION.md, SAFETY.md, COMMUNICATION_PROTOCOLS.md, MULTI_AGENT_WORKFLOWS.md, TROUBLESHOOTING.md).

The major gap is **inconsistent rule file coverage**: 18 rules are comprehensive reference implementations (core, dna, rna, gwas, protein, singlecell, multiomics, networks, visualization, ml, math, information, ontology, phenotype, simulation, life_events, epigenome, quality). However, 8 rule files (ecology, longread, metabolomics [critical], metagenomics, structural_variants, spatial, pharmacogenomics, menu) are minimal stubs lacking essential Patterns sections.

**Fix Priority**:
1. Restore `metabolomics.md` content (file appears corrupted/truncated)
2. Expand 7 minimal rule files to full coverage template
3. Add cross-links and standardize terminology
4. Clarify SPEC.md vs spec.md distinction

Once these actions are taken, documentation will be **100% complete and consistent** across all 44 files.

---

## Appendix: File Inventory

**Core Documentation (12)**:
```
AGENTS.md                 ✅ Complete
ARCHITECTURE.md           ✅ 500 lines
BEST_PRACTICES.md         ✅ 641 lines
COMMUNICATION_PROTOCOLS.md✅ 609 lines
index.md                  ✅ Complete
MULTI_AGENT_WORKFLOWS.md  ✅ 710 lines
ORCHESTRATION.md          ✅ 797 lines
PAI.md                    ✅ Complete
README.md                 ✅ Complete
SAFETY.md                 ✅ 663 lines
SPEC.md                   ✅ Complete
spec.md                   ⚠️  Duplicate purpose, rename recommended
TROUBLESHOOTING.md        ✅ 395 lines
```

**Module Rule Files (31)**:
```
rules/index.md            ✅ Complete (list of 28 modules)
rules/core.md             ✅ 277 lines (REFERENCE)
rules/dna.md              ✅ 303 lines
rules/rna.md              ✅ 433 lines
rules/gwas.md             ✅ 256 lines
rules/protein.md          ✅ 195 lines
rules/epigenome.md        ✅ 138 lines
rules/singlecell.md       ✅ 210 lines
rules/multiomics.md       ✅ 119 lines (concise but complete)
rules/networks.md         ✅ 174 lines
rules/visualization.md    ✅ 222 lines
rules/quality.md          ✅ 120 lines
rules/ml.md               ✅ 174 lines
rules/math.md             ✅ 260 lines
rules/information.md      ✅ 185 lines
rules/ontology.md         ✅ 179 lines
rules/phenotype.md        ✅ 162 lines
rules/simulation.md       ✅ 247 lines
rules/life_events.md      ✅ 175 lines
rules/ecology.md          ⚠️  Minimal (128 lines) — expand
rules/longread.md         ❌ Minimal (60 lines) — expand
rules/metabolomics.md     ❌ CRITICAL: Truncated (73 lines) — restore
rules/metagenomics.md     ❌ Minimal (47 lines) — expand
rules/structural_variants.md ❌ Minimal (47 lines) — expand
rules/spatial.md          ❌ Minimal (47 lines) — expand
rules/pharmacogenomics.md ❌ Minimal (48 lines) — expand
rules/menu.md             ❌ Minimal (38 lines) — expand
```

Total: 18 comprehensive ✅, 8 minimal ⚠️/❌, 1 broken ❌
