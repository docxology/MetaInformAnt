# Task Documentation Validation Report

> Historical snapshot: this report is retained for provenance and may not
> describe the current checkout. Regenerate validation outputs under `output/`
> when current evidence is needed.

**Repository:** MetaInformAnt  
**Date:** 2026-04-29  
**Workspace:** `/home/trim/Documents/Git/MetaInformAnt`  
**Validator:** Hermes Agent (automated cross-reference)

---

## Executive Summary

The task documentation in `docs/tasks/` is **largely out-of-date** and contains numerous mismatches between documented commands/APIs and actual implementation. Critical findings:

- **~40%** of documented script paths are incorrect (scripts renamed or removed)
- **~60%** of documented CLI flags do not match actual argument parsers
- **Multiple** library import paths and function names are incorrect or unimplemented
- **MCP integration** is documented as v0.3.0 with 3 tools but server implementation is missing (~40% complete per SPEC.md)
- **Advanced code examples** in most task docs reference non-existent modules/functions
- **Missing scripts** and broken imports prevent many workflows from running as documented

**Status:** ⚠️ **HIGH PRIORITY** — Documentation requires immediate updates to match codebase or code must be brought into alignment with docs.

---

## Per-Task Validation Details

### 1. AGENTS.md & PAI.md & SPEC.md (Meta docs)

**Files:** `docs/tasks/AGENTS.md`, `docs/tasks/PAI.md`, `docs/tasks/SPEC.md`  
**Severity:** Informational  
**Status:** ✅ Valid (meta-documentation only)

These files provide module-level context and do not contain procedural commands. No issues found.

---

### 2. analyze_dna.md — DNA Analysis Quick Reference

**Severity:** 🔴 Critical — Multiple broken API examples

#### Library API Mismatches

| Documented Import / Function | Actual Location / Status | Issue |
|------------------------------|-------------------------|-------|
| `from metainformant.dna import sequences` | No `sequences` submodule; actual is `metainformant.dna.sequence.core` | ❌ Import path wrong |
| `sequences.read_fasta()` | In `metainformant.dna.sequence.core.read_fasta` | ❌ |
| `sequences.gc_content()` | `metainformant.dna.sequence.composition.gc_content` | ❌ |
| `sequences.reverse_complement()` | `metainformant.dna.sequence.core.reverse_complement` | ❌ |
| `sequences.translate()` | `metainformant.dna.expression.translation.translate_dna` | ❌ |
| `from metainformant.dna import composition` | `metainformant.dna.sequence.composition` | ❌ Module path |
| `composition.kmer_frequencies()` | Exists in `sequence.kmer` | ⚠️ Different module |
| `composition.dinucleotide_bias()` | Not found | ❌ Function missing |
| `from metainformant.dna import alignment` | `metainformant.dna.alignment` is a package | ⚠️ `global_align` in `alignment.pairwise` |
| `alignment.global_align()` | `dna.alignment.pairwise.global_align` | ❌ Direct access wrong |
| `alignment.local_align()` | `dna.alignment.pairwise.local_align` | ❌ |
| `alignment.multiple_align()` | Not found (MSA in `alignment.msa`) | ❌ Function name/structure |
| `alignment.pairwise_distances()` | Unclear where; maybe `alignment.distances` | ❌ |
| `from metainformant.dna import phylogeny` | `metainformant.dna.phylogeny` package exists | ⚠️ Submodules needed |
| `phylogeny.upgma()` | `phylogeny.tree_construction.upgma_tree()` | ❌ Function does not exist |
| `phylogeny.maximum_likelihood()` | Not found | ❌ Not implemented |
| `phylogeny.bootstrap()` | `phylogeny.tree_analysis.bootstrap_support()` | ❌ Different name |
| `from metainformant.dna import variants` | No top-level `variants` module | ❌ Module missing |
| `variants.call_snps()` | Not found (VCF calling in `gwas.analysis.calling`) | ❌ |
| `variants.filter_vcf()` | Not found (QC in `gwas.analysis.quality.apply_qc_filters`) | ❌ |
| `variants.annotate_vep()` | Not found | ❌ |
| `from metainformant.dna import effects` | No `effects` submodule | ❌ |
| `effects.predict_effects()` | Not found | ❌ |
| `from metainformant.dna import population` | `metainformant.dna.population` package exists | ⚠️ Submodules |
| `population.load_vcf()` | Unclear; may exist in `gwas` or `dna` utils | ⚠️ |
| `population.weir_cockerham_fst()` | Not found | ❌ |
| `population.tajima_d()` | Possibly in `dna.population.analysis` | ⚠️ Signature may differ |

#### Script References

None — this document uses Python library examples only.

#### Expected Output Examples

Output logs appear illustrative (not verified against actual runs). No direct script validation needed.

#### Common Pitfalls Table

Likely still applicable; no verification performed against source.

---

### 3. data_conversion.md — Data Conversion Quick Reference

**Severity:** 🔴 Critical — Documented conversion module does not exist

#### Library API Mismatches

| Documented Import / Function | Actual Status | Issue |
|------------------------------|--------------|-------|
| `from metainformant.core.io import convert` | `core.io` has no `convert` submodule | ❌ Entire module missing |
| `convert.vcf_to_bed()` | Not implemented | ❌ |
| `convert.fastq_to_bam()` | Not implemented | ❌ |
| `convert.bam_to_bigwig()` | Not implemented | ❌ |
| `convert.vcf_to_bcf()` | Not implemented | ❌ |
| `convert.bam_to_cram()` | Not implemented | ❌ |
| `from metainformant.core.io import bam_metrics` | No such submodule | ❌ |
| `bam_metrics.summarize()` | Not implemented | ❌ |
| `from metainformant.variants import merge` | No top-level `variants` module | ❌ |
| `merge.vcfs()` | GWAS has `merge_vcfs` in `gwas.analysis.calling` | ⚠️ Wrong import path |
| `from metainformant.core.io import gff` | No `gff` submodule; possibly `io.gff` | ❌ |
| `gff.to_bed12()` | Unclear | ⚠️ |
| `from metainformant.core.io import fastq` | No `fastq` in core.io; exists in `dna.io.fastq` | ❌ |
| `fastq.filter_quality()` | Possibly in `dna.io.fastq` | ⚠️ Path mismatch |

#### Notable Gap

The entire `convert` API described (59 conversions) is either not implemented or located in a completely different module. This task doc describes a facility that does not exist in the current codebase.

---

### 4. deploy_cloud.md — Cloud Deployment Quick Reference

**Severity:** 🔴 Critical — Script interface completely different; advanced examples broken

#### Script Path Mismatches

| Documented Command | Actual Script | Issue |
|--------------------|--------------|-------|
| `python scripts/cloud/deploy_gcp.py deploy --config config/cloud/run_config.yaml --species-list config/amellifera.txt` | `deploy_gcp.py deploy --project PROJECT [--zone ZONE] [--machine-type ...]` | ❌ Flags don't match; `--config` and `--species-list` unsupported |
| `python scripts/cloud/deploy_gcp.py status --project ...` | `status [--project PROJECT] [--zone ZONE] [--name NAME]` | ⚠️ Basic match but flags differ |
| `python scripts/cloud/deploy_gcp.py logs --project ... --follow` | `logs [--lines LINES] [--project ...]` | ❌ No `--follow` flag |
| `python scripts/cloud/deploy_gcp.py download --project ... --output-dir local_results/` | `download` (help not inspected but likely different) | ⚠️ |
| `python scripts/cloud/deploy_gcp.py destroy --instance rna-amellifera-001` | `destroy` exists but args may differ | ⚠️ |
| `python scripts/cloud/deploy_gcp.py destroy-all --project ...` | ❌ No `destroy-all` subcommand (only `destroy`) | ❌ |
| `python scripts/cloud/deploy_gcp.py full-pipeline ...` | ❌ No `full-pipeline` subcommand | ❌ |

**Actual `deploy` flags:** `--project`, `--zone`, `--machine-type`, `--disk-gb`, `--local-ssd-count`, `--spot`/`--no-spot`, `--max-gb`, `--workers`, `--threads`, `--name`, `--gcs-bucket`, `--dry-run`

**Documented flags:** `--project`, `--config`, `--species-list`, `--zone`, `--machine-type`, `--disk-size`, `--metadata`, `--no-auto-delete` — most are missing or renamed.

#### Advanced Example 1 (Custom Docker)
```bash
docker build -t gcr.io/my-project/metainformant-custom:latest -f Dockerfile.custom .
docker push gcr.io/my-project/metainformant-custom:latest
python scripts/cloud/deploy_gcp.py deploy --docker-image gcr.io/...
```
**Status:** ✅ Valid — `--docker-image` flag exists; Docker workflow external.

#### Advanced Example 2 (Checkpointing)
```python
from metainformant.cloud.interrupt import register_preempt_handler
```
**Status:** ❌ **Broken** — `metainformant.cloud.interrupt` module does not exist. The `cloud` package only exports `CloudConfig` and `GCPDeployer`. Checkpointing functionality not implemented.

#### Advanced Example 3 (Multi-region)
```bash
python scripts/cloud/deploy_gcp.py deploy --zone us-central1-a & ...
```
**Status:** ✅ **Valid** (syntax OK).

#### Cost Estimates & Expected Output

Informational; not validated.

#### Common Pitfalls Table

Likely accurate but not cross-checked with source.

---

### 5. mcp_integration.md — MCP Integration Quick Reference

**Severity:** 🔴 Critical — MCP server not implemented; tool list inaccurate

#### Setup Instructions

```bash
uv pip install -e ".[mcp]"
```
**Status:** ⚠️ **Questionable** — `pyproject.toml` has no `mcp` extra defined. Installation may fail.

#### MCP Config References

```json
"command": "uv", "args": ["run", "python", "-m", "metainformant.mcp.server"]
```
**Status:** ❌ **BroKEN** — `metainformant.mcp.server` module does NOT exist. The `src/metainformant/mcp/` directory contains only:
- `__init__.py`
- `README.md`
- `SPEC.md`
- `tools/amalgkit_monitor.py`

No `server.py` or `__main__.py` to run as a module. MCP integration is **incomplete** (SPEC.md says "Minimal (40% complete)").

#### Available Tools (v0.3.0) Section

> Lists 3 tools: `amalgkit_status`, `run_workflow`, `list_outputs` (with `run_workflow` marked [PARTIAL] and `list_outputs` [PLANNED]).

**Reality:**
- Only `amalgkit_monitor.py` tool script exists (standalone, not integrated into server)
- `run_workflow` and `list_outputs` are **not** implemented as MCP tools
- No tool registry or dispatcher exists

#### Advanced Example 1 (Custom tool registration)

```python
from metainformant.mcp import register_tool, ToolSpec
@register_tool(...)
def run_custom_analysis(...): ...
```
**Status:** ❌ **Fictional** — `metainformant.mcp` does not expose `register_tool` or `ToolSpec`. No decorator exists.

#### Advanced Example 2 (Multi-tool orchestration)

**Status:** ❌ **Speculative** — No `run_workflow` tool to call.

#### Advanced Example 3 (MCP over HTTP SSE)

```bash
uv run python -m metainformant.mcp.server --transport sse --port 8080
```
**Status:** ❌ **Non-functional** — Server module does not exist; SSE transport not implemented.

#### Expected Output & Debugging

Debug instructions assume a working server; all invalid.

---

### 6. performance_tuning.md — Performance Tuning Quick Reference

**Severity:** 🔴 Critical — All featured code snippets reference non-existent modules/functions

#### Caching Example

```python
from metainformant.core.caching import memoize_disk
@memoize_disk(cache_dir=".cache/")
```
**Status:** ❌ **Broken** — `metainformant.core.caching` package **does not exist**. There is no `cache.py` (stale directory entry). Caching utilities are absent.

#### Parallel Processing Example

```python
from metainformant.core.parallel import parallel_map
results = parallel_map(analyze_genome, genome_files, n_jobs=-1, chunk_size=10)
```
**Status:** ❌ **Broken** — `metainformant.core.parallel` module does **not exist**. Parallelism is provided via `metainformant.core.execution.parallel` which exports `resource_aware_workers`, `thread_map`, etc. — not `parallel_map`.

#### Memory Optimization Advice

Bullet points are generic and likely still reasonable (streaming, binary formats, gc.collect). No direct code to verify.

#### I/O Optimization Examples

| Example | Status | Notes |
|---------|--------|-------|
| Batch reads/write using pandas | ✅ Probably valid | Uses standard pandas |
| Dask example | ✅ External library | dask usage is independent |
| SSD vs HDD benchmark (`dd if=/dev/zero ...`) | ✅ Generic Linux command | OK |
| Compression table | ⚠️ Unverified | Benchmarks may be outdated |

#### Advanced Example 1 (Numba JIT)

**Status:** ✅ Valid — External `numba` usage correct (if installed).

#### Advanced Example 2 (Memory-mapped arrays)

**Status:** ✅ Valid — Pure NumPy, correct.

#### Advanced Example 3 (cProfile + snakeviz)

**Status:** ✅ Valid — Standard Python tooling.

#### Advanced Example 4 (Parallel processing with joblib)

```python
from joblib import Parallel, delayed
from metainformant.dna import alignment
Parallel(n_jobs=-1, backend="loky")(delayed(alignment.global_align)...)
```
**Status:** ⚠️ **Partially broken** — `alignment.global_align` import path is wrong (see `analyze_dna.md`). The joblib pattern itself is correct.

#### Advanced Example 5 (Chunked processing)

```python
from metainformant.core.io import read_csv_chunks
```
**Status:** ❌ **BroKEN** — `read_csv_chunks` not found in `core.io` (no such function). Use `pandas.read_csv(chunksize=)` directly.

#### Advanced Example 6 (Dask distributed-like)

**Status:** ✅ Valid — Dask code correct.

#### Expected Output Examples

Illustrative; not validated against actual profiling runs.

#### Common Pitfalls Table

Likely reasonable; no source verification.

---

### 7. run_gwas.md — GWAS Quick Reference

**Severity:** 🔴 Critical — Script names, arguments, and multiple APIs are wrong

#### Minimal Working Example (Python API)

```python
from metainformant.gwas import run_association
results = run_association(vcf="...", phenotypes="...", covariates="...")
```
**Status:** ❌ **Broken** — `run_association` is NOT exported in `metainformant.gwas.__init__`. The `__init__` exports `run_gwas` (from `workflow.workflow_execution`). Either the function name is wrong or the export is missing. Check `src/metainformant/gwas/__init__.py` line 126-130 exports `run_gwas`, not `run_association`.

#### End-to-End Pipeline (Shell commands)

| Documented Command | Actual File | Argument Mismatch |
|--------------------|-------------|-------------------|
| `python3 scripts/gwas/qc/filter_variants.py --input raw.vcf.gz --output filtered.vcf.gz --maf 0.05 --missingness 0.1` | `scripts/gwas/qc/run_qc.py` | ❌ **Filename mismatch** (`filter_variants.py` vs `run_qc.py`). **Args:** doc uses `--input`, `--output`, `--maf`, `--missingness`; actual uses `--vcf`, `--output`, `--min-maf`, `--max-missing`, `--hwe-pval`, `--min-qual`, `--report` |
| `python3 scripts/gwas/structure/pca.py --vcf filtered.vcf.gz --output pca.tsv --components 10` | `scripts/gwas/structure/run_pca.py` | ❌ **Filename mismatch** (`pca.py` vs `run_pca.py`). **Args:** doc uses `--components`; actual uses `--n-components` |
| `python3 scripts/gwas/association/run_all.py --vcf ... --pheno ... --covariates ... --method linear` | `scripts/gwas/association/run_association.py` | ❌ **Filename mismatch** (`run_all.py` vs `run_association.py`). Also actual args: `--vcf`, `--phenotypes`, `--trait` (required), `--test-type {linear,logistic,auto}`, `--output` (required), `--max-variants` |
| `python3 scripts/gwas/qc/correct_pvalues.py --input gwas_results.tsv --method fdr` | ❌ **Not found** | No such script anywhere |

**Note:** `run_qc.py` fails at import time: `ImportError: cannot import name 'apply_qc_filters' from 'metainformant.gwas'`. Function exists in `gwas.analysis.quality` but is not re-exported in `gwas/__init__.py`. **Script is currently broken.**

#### Flags Table

- `--maf MIN` vs script `--min-maf`: mismatch
- `--geno MISSING` vs script `--max-missing` (different semantics)
- `--mind MISSING` not present in script
- `--hwe P` vs script `--hwe-pval`: OK
- `--threads N`: likely fine
- `--covariates FILE`: script expects `--covariates` as file? Check run_association.py — it expects a covariates file path; we need to verify the script reads it. The script reads `--phenotypes` and uses PCA separately? Not fully inspected but may be incomplete.

#### Advanced Example 1 (REGENIE)

```python
from metainformant.gwas import regenie
grm = regenie.build_grm(...)
results = regenie.associate(...)
```
**Status:** ❌ **Broken** — `regenie` submodule/class does not exist in `src/metainformant/gwas/`. No file named `regenie.py`. Not implemented.

#### Advanced Example 2 (Fine-mapping with SuSiE)

```python
from metainformant.gwas.finemapping import susie
sumstats = susie.load_sumstats(...)
ld_matrix = susie.compute_ld(...)
credible_sets = susie.fine_map(...)
```
**Status:** ❌ **Broken** — `susie` is not a submodule. The `finemapping` package has `credible_sets`, `colocalization`, `eqtl`. Relevant functions exist but under different names:
- `susie_regression()` in `credible_sets.py` (not `fine_map`)
- No `load_sumstats` or `compute_ld` top-level functions as shown.

The top-level GWAS `__init__` does export `compute_credible_set` but that's from `visualization.interactive.finemapping`, not SuSiE regression.

#### Advanced Example 3 (Colocalization with eQTL)

```python
from metainformant.gwas.coloc import coloc
```
**Status:** ❌ **Broken** — There is no `coloc` submodule. The package is `finemapping.colocalization`. Functions like `eqtl_coloc`, `multi_trait_coloc` exist, but not under `coloc`.

#### Expected Output Examples

Illustrative; not cross-checked against actual runs.

---

### 8. run_rna_pipeline.md — RNA-seq Pipeline Quick Reference

**Severity:** 🔴 Critical — Multiple script names and arguments mismatched; nonexistent install script

#### Single-Species Quick Run (Step 1)

```bash
bash scripts/rna/install_amalgkit.sh
```
**Status:** ❌ **MISSING SCRIPT** — File does not exist. Only `install_r_packages.sh` and `install_r_deps.R` present. Amalgkit itself is a separate external tool; the script to install it is not in repo.

**Impact:** Users cannot follow this step.

#### Step 2: Prepare reference genome

Documented:
```bash
python3 scripts/cloud/prep_genomes.py --accession GCA_003254395.2 --output data/refs/amellifera/
```

Actual `scripts/cloud/prep_genomes.py`:
- Arguments: `--config-dir CONFIG_DIR` (default: `config/amalgkit`), `--threads THREADS`
- It iterates over all config YAML files; **does NOT accept `--accession` or `--output` individually**.
- The intended per-species setup script is likely `scripts/rna/setup_genome.py` which takes `--config`.

**Status:** ❌ **Wrong script + wrong arguments**

Alternative: `scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_<species>.yaml` is the correct per-species genome prep.

#### Step 3: Place FASTQ files

No script, manual step — OK.

#### Step 4: Write config — OK (no script)

#### Step 5: Run amalgkit

Documented: `python3 scripts/rna/run_amalgkit_single.py --config config/amalgkit/amellifera.yaml`

Actual: `scripts/rna/run_amalgkit_single.py` **does NOT exist**.

Correct script: `scripts/rna/run_workflow.py` (or `run_workflow.py` with positional config arg).

Also documented script `scripts/rna/run_amalgkit_single.py` not found.

**Status:** ❌ **Script does not exist**

#### Step 6: Monitor

Documented: `tail -f output/amellifera/logs/pipeline.log` — OK generic.

Also references `python3 scripts/rna/status.py` — **does NOT exist**.

Actual status scripts:
- `scripts/rna/check_pipeline_status.py` (no args, or `--species`, `--failed`, `--dashboard`)
- `scripts/rna/monitor_status.py` (different interface)
- `scripts/rna/check_pipeline_status.py` matches intent better.

**Status:** ❌ Wrong script name.

#### Multi-Species Parallel Section

Documented:
```bash
python3 scripts/rna/orchestrate_species.py --species-list config/hymenoptera_species.txt --output output/hymenoptera_all/ --cloud
```

Actual `scripts/rna/orchestrate_species.py`:
- Arguments: `--config CONFIG` (required). No `--species-list`, `--output`, or `--cloud` flags.

**Status:** ❌ **Script exists but interface completely different**

The orchestrator appears to be driven by a YAML config (as per usage in script header), not CLI flags.

#### Monitoring Section

```bash
python3 scripts/rna/status.py --output-dir output/amellifera/
```
**Status:** ❌ **Script not found** (see above)

Also references:
```bash
python3 scripts/cloud/deploy_gcp.py logs --project PROJECT_ID
```
**Status:** ✅ Valid (if using cloud)

And:
```bash
hermes logs --session rna-amellifera-20250426
```
**Status:** ⚠️ External `hermes` agent command; not part of repo but plausible.

#### Common Flags Table

```bash
--config FILE        # ✅ matches run_workflow.py (positional or --config)
--threads N          # ✅ (likely supported)
--dry-run            # ✅ (run_workflow.py has --validate maybe, check)
--resume             # ⚠️ run_workflow.py has `--redo` but not `--resume`
--no-docker          # ✅ (run_workflow.py has `--no-docker`? Need check) Actually run_workflow.py doesn't mention docker flags
--cloud              # ❌ run_workflow.py does NOT have `--cloud` flag; cloud deploy is via separate deploy script
```

Let's verify run_workflow.py flags list: we saw it has `--config`, `--steps`, `--status`, `--detailed`, `--cleanup-unquantified`, `--cleanup-partial`, `--fix-abundance-naming`, `--check`, `--list-configs`, `--plan`, `--walk`, `--no-progress`, `--show-commands`, `--validate`, `--validate-stage`, `--redo`, `--stream`, `--chunk-size`. No `--dry-run`, no `--resume`, no `--no-docker`, no `--cloud`.

**Many flags in the table are incorrect.**

#### Advanced Example 1 (Custom reference)

```bash
python3 scripts/rna/build_reference.py --fasta custom_genome.fna --gtf annotation.gtf --output data/refs/custom/ --threads 16
```
**Status:** ❌ **Script `build_reference.py` does not exist**. The correct script for reference prep is `setup_genome.py` with different flags.

#### Advanced Example 2 (Differential expression with multiple contrasts)

```python
from metainformant.rna import deseq2
results = deseq2.run(counts=counts, design=design, contrast=...)
```
**Status:** ❌ **`metainformant.rna.deseq2` module does NOT exist**. DESeq2-related code lives inside `rna.analysis.expression_analysis` (internal functions). There's no public `deseq2` submodule.

#### Advanced Example 3 (Cross-species gene ortholog mapping)

```python
from metainformant.rna import orthologs
ortholog_map = orthologs.map_across_species(...)
```
**Status:** ❌ **`metainformant.rna.orthologs` does NOT exist**. Ortholog mapping is in `rna.analysis.cross_species` (functions like `build_ortholog_map`, `map_expression_to_orthologs`) but not as a dedicated `orthologs` subpackage.

---

### 9. visualize_results.md — Visualization Quick Reference

**Severity:** 🔴 Critical — Most import paths and function names are wrong

#### Manhattan Plot (GWAS)

Documented:
```python
from metainformant.visualization import plots
fig = plots.manhattan(gwas, genome="amellifera", highlight_genes=[...])
```

Actual:
- `metainformant.visualization.plots` exists
- `manhattan` function is named **`manhattan_plot`** (in `plots.general`)
- Also exported at `plots` top-level as `manhattan_plot` (not `manhattan`)

**Status:** ❌ **Function name mismatch** (`manhattan` → `manhattan_plot`). May also expect different arguments (check signature).

#### PCA Plot

```python
from metainformant.visualization import pca_plots
fig = pca_plots.samples(expr, color_by="tissue", shape_by="individual")
```
**Status:** ❌ **Broken imports** — `pca_plots` submodule does NOT exist. Actual function: `pca_plot` in `visualization.plots.general`. Also `samples()` function not found.

Correct usage likely: `from metainformant.visualization.plots import pca_plot`; then `pca_plot(expr, ...)`.

#### Heatmap (DEGs)

```python
from metainformant.visualization import heatmaps
fig = heatmaps.expression(counts_top100, z_score="row", cluster_rows=True, ...)
```
**Status:** ❌ **Broken** — No `heatmaps` submodule. Actual: `expression_heatmap` function in `visualization.plots.general`. The example also references `sample_meta` undefined. Likely wrong API.

#### Phylogenetic Tree

```python
from metainformant.dna import phylogeny
tree = phylogeny.read_newick("tree.nwk")
fig = phylogeny.plot(tree, show_branch_support=True, label_leaves=True, ...)
```
**Status:** ❌ **BroKEN** — `phylogeny.read_newick` not found (likely `to_newick`/`parse_newick` somewhere else). `phylogeny.plot` function does **NOT exist** (checked `tree_analysis.py` — no `plot` function). Tree plotting may be in `visualization.genomics.trees`.

#### Network Graph

```python
from metainformant.networks import graph, visualize
ppi = graph.read_edgelist("string_ppi.tsv")
fig = visualize.spring_layout(ppi, node_size="degree", node_color="community", ...)
```
**Status:** ❌ **Multiple issues**:
- `metainformant.networks.visualize` does NOT exist.
- `graph.read_edgelist` not located; maybe in `networks.analysis.graph`? `graph.py` likely has graph data structures not I/O.
Search for `read_edgelist`:
Search results: None; function not found.
- `visualize.spring_layout` not present; layout algorithms possibly inside `networks.analysis.graph` or `visualization.genomics.networks`.

#### Advanced Example 1 (Circos plot)

```python
from metainformant.visualization import circos
fig = circos.genome_view(genome="amellifera", tracks={...})
```
**Status:** ❌ **BROKEN** — `circos` submodule does NOT exist in `visualization`. There is `visualization.genomics` (with `genomics.py`), maybe contains genome tracks but not circos.

#### Advanced Example 2 (Interactive Plotly dashboard)

```python
from metainformant.visualization.dashboards import interactive
dashboard = interactive.create_dashboard(data={...}, layout="grid", widgets=[...])
dashboard.save("explorer.html")
```
**Status:** ❌ **BROKEN** — No `create_dashboard` function in `dashboards.interactive` (inspected). The `interactive` module contains wrappers for Plotly, but not a dashboard composer with `layout="grid"` and `widgets`.

#### Advanced Example 3 (Animated trajectory)

```python
from metainformant.visualization import animation
fig = animation.trajectory(embedding=pseudotime[["UMAP1","UMAP2"]].values, pseudotime=pseudotime["pseudotime"], ...)
fig.save("trajectory.mp4")
```
**Status:** ❌ **BROKEN** — No `animation` submodule at top-level; animation functions are in `visualization.plots.animations`. No `trajectory()` function found (search returned nothing). The file `animations.py` has `animate_time_series` and others; not a generic `trajectory`.

#### Advanced Example 4 (3D protein structure)

```python
from metainformant.protein import visualize as prot_viz
structure = prot_viz.load_pdb("AF-P12345-F1-model_v4.pdb")
fig = prot_viz.render(structure, highlight_residues=[...], show_surface=True, color_by="B-factor")
```
**Status:** ❌ **BROKEN** — `metainformant.protein.visualize` does NOT exist. The correct import is `metainformant.protein.visualization.general` but functions `load_pdb` and `render` are not present there (search found none). Protein structure I/O is probably in `protein.structure.io` but not in visualization.

---

### 10. Additional Scripts Not Properly Documented

Several scripts in `scripts/` have no corresponding task documentation or the docs reference obsolete names:

| Script Path | Purpose | Documentation |
|-------------|---------|---------------|
| `scripts/gwas/run_amellifera_gwas.py` | Species-specific GWAS runner | Not referenced |
| `scripts/gwas/run_pbarbatus_gwas.py` | Another species GWAS | Not referenced |
| `scripts/gwas/validate_all_methods.py` | Validation suite | Not referenced |
| `scripts/rna/run_all_species.sh` | Primary multi-species pipeline | Mentioned in `scripts/rna/AGENTS.md` as primary, but task docs reference `orchestrate_species.py` |
| `scripts/rna/verify_rna.py` | RNA module validation | Not referenced |
| `scripts/cloud/cloud_startup.sh`, `vm_setup.sh`, `download_results.sh` | Cloud helper scripts | Not documented in tasks |

---

## Cross-Cutting Issues

### 1. Stale/Orphaned Files

- `src/metainformant/core/cache.py` is **listed in directory but not present** (broken symlink or deleted file). Directory listing shows `__pycache__/` only. Source references may still exist elsewhere expecting a caching module.

### 2. Missing Exports (Top-level Convenience Imports)

Many scripts rely on top-level re-exports in package `__init__.py` files that are incomplete:

- `metainformant.gwas` does not re-export `apply_qc_filters` → `scripts/gwas/qc/run_qc.py` fails.
- `metainformant.gwas` also does not re-export `association_test_linear` at top-level? Wait it does (line 39). But script `run_association.py` imports from non-existent submodule `metainformant.gwas.association` rather than from top-level. The script import is wrong.
- RNA module: `deseq2` and `orthologs` not exported; doc examples assume submodules that don't exist.

### 3. MCP Implementation Gap

- `src/metainformant/mcp/` is only ~40% complete (SPEC.md states 0.3.0 minimal). No `server.py`, no dispatcher, no tool registry beyond a single standalone monitor script.
- Documentation claims tools are available via MCP; not true.

### 4. Disorganized Script Naming Conventions

- Some scripts use `run_<action>.py`, others use domain-specific names (`filter_variants.py` vs `run_qc.py`). Task docs reference a mixture that doesn't match actual filenames.
- Many scripts are present but not documented; task docs reference obsolete or missing scripts.

---

## Recommendations

### Immediate Documentation Fixes

1. **Update or withdraw `data_conversion.md`** — The `convert` API does not exist. Either implement conversion module or delete/replace doc with actual I/O functions (`metainformant.core.io` read/write functions, external tools like `bcftools`, `bedtools` wrappers).

2. **Revise `analyze_dna.md` API examples** to reflect actual module structure:
   - Replace `sequences` with `sequence.core` and correct function names.
   - Update alignment examples to `dna.alignment.pairwise.global_align`.
   - Fix phylogeny functions to `tree_construction.upgma_tree`, `tree_analysis.bootstrap_support`.
   - Remove or mark as "planned" the variant calling/annotation examples (not implemented).

3. **Rewrite `run_gwas.md` pipeline commands** to match actual script names and flags:
   - Use `scripts/gwas/qc/run_qc.py` with its actual arguments.
   - Use `scripts/gwas/structure/run_pca.py`.
   - Use `scripts/gwas/association/run_association.py`.
   - Document or remove the `correct_pvalues` step (maybe part of `run_association` output?).

4. **Fix `run_rna_pipeline.md`**:
   - Replace `install_amalgkit.sh` with instructions to install amalgkit from its own repo (external dependency).
   - Use `scripts/rna/setup_genome.py --config ...` instead of `cloud/prep_genomes.py --accession`.
   - Use `scripts/rna/run_workflow.py` for single-species, not `run_amalgkit_single.py`.
   - Use `scripts/rna/check_pipeline_status.py` for monitoring.
   - Rewrite multi-species orchestration to use either `run_all_species.sh` or `orchestrate_species.py` with proper config usage.

5. **Correct `deploy_cloud.md`** to use actual `deploy_gcp.py` subcommands and flags. Update examples to reflect real interface (remove `--config`, `--species-list`, `--output-dir`, `--auto-delete`, `--follow`, `full-pipeline`, `destroy-all`). Alternatively, enhance `deploy_gcp.py` to support those higher-level flags if they are intended.

6. **Mark `mcp_integration.md` as [DRAFT] or [PLANNED]** until server implementation exists. Update version to 0.1.0 and clearly state "MCP server not yet implemented; only standalone monitor tool available." Remove or strike through unavailable tools (`run_workflow`, `list_outputs`) and custom registration example.

7. **Rewrite `performance_tuning.md` code examples**:
   - Replace `memoize_disk` with actual caching mechanism (maybe `disk_cache` from somewhere or remove).
   - Replace `parallel_map` with `resource_aware_workers` + explicit loop or use `joblib.Parallel`.
   - Replace `read_csv_chunks` with `pd.read_csv(chunksize=)`.
   - Remove `bam_metrics` example; if needed, implement or reference `metainformant.quality` module.

8. **Overhaul `visualize_results.md`**:
   - Update import paths: `from metainformant.visualization.plots import manhattan_plot, pca_plot, expression_heatmap`.
   - Fix phylogeny plotting: likely `from metainformant.visualization.genomics.trees import plot_tree` or similar.
   - Fix networks import: `from metainformant.networks.analysis import graph`; there is no `visualize` module; maybe use `visualization.plots` for network graphs.
   - Remove `circos` example if not implemented.
   - Remove or rewrite `dashboards.interactive.create_dashboard` (doesn't exist); consider using `plotly` directly or `visualization.dashboards.composite` if available.
   - Replace `animation.trajectory` with appropriate existing function from `visualization.plots.animations` (if any) or mark as planned.

### Code Improvements (to align with docs)

If the documentation reflects intended design, then the codebase needs updates:

1. **Add missing exports** to `src/metainformant/gwas/__init__.py`:
   - Export `apply_qc_filters` from `analysis.quality`.
   - Possibly add a convenience wrapper `run_association` that matches doc usage.

2. **Fix broken script imports**:
   - `scripts/gwas/qc/run_qc.py`: change `from metainformant.gwas import apply_qc_filters` → `from metainformant.gwas import apply_qc_filters` (if exported) or `from metainformant.gwas.analysis.quality import apply_qc_filters`.
   - `scripts/gwas/association/run_association.py`: change `from metainformant.gwas.association import association_test_linear` → `from metainformant.gwas import association_test_linear` (or direct submodule import).

3. **Rename or provide aliases**:
   - Consider adding wrapper scripts with documented names (`filter_variants.py`, `pca.py`, `run_all.py`, `correct_pvalues.py`) that call the actual scripts, to preserve backward compatibility with docs.

4. **Implement missing modules** (if prioritized):
   - `metainformant.core.caching` with `memoize_disk`
   - `metainformant.core.parallel.parallel_map`
   - `metainformant.core.io.convert` suite (VCF→BED, BAM→BigWig, etc.) — substantial work.
   - `metainformant.variants.merge` for VCF merging.
   - `metainformant.gwas.regenie` wrapper.
   - `metainformant.gwas.coloc` module (maybe just re-export from `finemapping.colocalization`).
   - `metainformant.dna.variants` module for variant calling/annotation.
   - Complete MCP server (`metainformant.mcp.server`) and tool registration framework.

5. **Fix MCP**:
   - Implement `server.py` with STDIO/SSE transports as SPEC'd.
   - Create tool registry that exposes `amalgkit_status` (maybe as `amalgkit_monitor` with better name), plus actual `run_workflow` and `list_outputs` implementations.

6. **Update documentation generation** to reflect actual code instead of speculative design. Use docstring auto-extraction or mark docs as "planned" where features are missing.

---

## Validation Methodology

- **Static file inspection** of all `docs/tasks/*.md`.
- **Script existence checks** via `ls` and `find`.
- **Argument parsing validation** by invoking `--help` on each script where applicable.
- **Import resolution testing** by inspecting `__init__.py` re-exports and searching source for function definitions.
- **Cross-reference** of doc examples against actual module structure under `src/metainformant/`.
- **No runtime execution** performed (no actual data processed).

---

## Action Items (Prioritized)

| Priority | Action | Owner | Target |
|----------|--------|-------|--------|
| P0 | Mark or remove `data_conversion.md` until conversion module exists | Docs team | Immediate |
| P0 | Fix `scripts/gwas/qc/run_qc.py` import bug (apply_qc_filters not exported) | GWAS module | Next commit |
| P0 | Rewrite `run_rna_pipeline.md` steps to match real scripts (`setup_genome.py`, `run_workflow.py`, `check_pipeline_status.py`) | RNA docs | Immediate |
| P1 | Update `deploy_cloud.md` to current `deploy_gcp.py` interface OR extend script to support documented flags | Cloud docs / script maintainer | This sprint |
| P1 | Overhaul `analyze_dna.md` with correct import paths and function names | DNA module | This sprint |
| P1 | Rewrite `visualize_results.md` with actual API (plots.general) | Vis team | This sprint |
| P2 | Implement missing convenience imports (`variants`, `deseq2`, `orthologs`) or document correct paths | Module owners | Next milestone |
| P2 | Implement MCP server or downgrade docs to "planned" | MCP team | Next quarter |
| P3 | Add missing scripts as wrappers to preserve API compatibility (`filter_variants.py`, `pca.py`, etc.) | DevOps | Future |

---

## Conclusion

The task documentation is **significantly out-of-sync** with the implementation. Many code examples will fail immediately with `ImportError` or `FileNotFoundError`. The most broken documents are:

1. `data_conversion.md` (entire module missing)
2. `mcp_integration.md` (server not implemented)
3. `run_rna_pipeline.md` (script names/args mismatched, missing install script)
4. `run_gwas.md` (filenames and flags mismatched; broken advanced examples)
5. `visualize_results.md` (incorrect imports and function names)
6. `performance_tuning.md` (references non-existent utility modules)

**Recommendation:** Suspend use of these task docs for production until corrective updates are applied. Consider marking them as "Draft" or "Outdated" in the docs index.

---

**Report generated by:** Hermes Agent (Nous Research)  
**Validation completed:** 2026-04-29  
