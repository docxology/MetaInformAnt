# RNA Steps

Each AMALGKIT subcommand is exposed via a dedicated runner in `metainformant.rna.steps`.

Benefits:

- Step-specific hooks live in isolated modules
- Stable call surface: `run(params, work_dir=..., log_dir=..., check=...)`
- Orchestrator remains thin and composable

## Example

```python
from pathlib import Path
from metainformant.rna.steps import STEP_RUNNERS

params = {"threads": 4}
runner = STEP_RUNNERS["metadata"]
res = runner(params, work_dir=Path("output/amalgkit/run1"), log_dir=Path("output/amalgkit/run1/logs"))
print(res.returncode)
```

## Available step modules

- metadata
- integrate
- config
- select
- getfastq
- quant
- merge
- cstmm
- curate
- csca
- sanity

## RNA: Step Runners

Each `amalgkit` subcommand has a corresponding runner in `metainformant.rna.steps` for internal orchestration. All share the same signature:

```python
def run(params: Mapping[str, Any] | None = None, *, work_dir: str | Path | None = None, log_dir: str | Path | None = None, check: bool = False) -> CompletedProcess[str]
```

Available runners: `metadata`, `integrate`, `config`, `select`, `getfastq`, `quant`, `merge`, `cstmm`, `curate`, `csca`, `sanity`.

```python
from pathlib import Path
from metainformant.rna.steps import STEP_RUNNERS

runner = STEP_RUNNERS["quant"]
res = runner({"threads": 4}, work_dir=Path("./work"), log_dir=Path("./work/logs"), check=False)
```

## Quick Links

- **[Step Documentation](amalgkit/steps/README.md)** - All 11 step guides with function signatures
- **[API Reference](API.md#step-runner-functions)** - Step runner function documentation
- **[Function Index](amalgkit/FUNCTIONS.md)** - Quick function lookup
- **[Workflow Guide](workflow.md)** - Workflow planning and execution
- **[Main Index](README.md)** - RNA domain master index

Related: [workflow.md](workflow.md), [amalgkit/README.md](amalgkit/README.md)

## See Also

### Documentation
- **[Step Documentation](amalgkit/steps/README.md)** - Detailed step guides with function signatures
- **[API Reference](API.md#step-runner-functions)** - Complete function documentation
- **[Function Index](amalgkit/FUNCTIONS.md)** - Quick function lookup
- **[Workflow Guide](workflow.md)** - Workflow planning and execution
- **[Configuration Guide](CONFIGURATION.md)** - Configuration management
- **[Orchestration Guide](ORCHESTRATION/README.md)** - Orchestrator overview

## Step-by-step details

The sections below summarize goals, common parameters, and expected artifacts. Parameters are passed as a mapping and converted to CLI flags in `amalgkit`; see `build_cli_args` rules in `metainformant.rna.amalgkit`.

### metadata
- Purpose: discover SRA/ENA accessions and associated sample metadata for target taxa/species/tissues
- Common params: `taxon-id` (int), `species-list` (list[str]), `tissue` (list[str]), `threads` (int)
- Artifacts: metadata tables under `work_dir` (exact paths defined by amalgkit); logs under `work_dir/logs`

### integrate
- Purpose: harmonize metadata from different sources into a consistent table for selection
- Common params: inherits common threads/species; additional filters may be provided per amalgkit
- Artifacts: integrated metadata table

### config
- Purpose: materialize per-run configuration for downstream steps
- Notes: often precedes selection and acquisition; no heavy I/O

### select
- Purpose: filter samples/runs by metadata criteria (e.g., tissues, platforms)
- Common params: predicates like `tissue`, study filters, min length/quality as supported by amalgkit
- Artifacts: selection list/use-manifest for downloads and quantification

### getfastq
- Purpose: download FASTQ files for selected runs
- Common params: `out-dir` (Path; default layout uses `work_dir/fastq`), `threads`
- Artifacts: FASTQ files organized by run/sample; logs in `logs/`

### quant
- Purpose: quantify transcript/gene abundances (e.g., Salmon)
- Common params: `out-dir` (default `work_dir/quant`), index/reference options as required by amalgkit, `threads`
- Artifacts: per-sample quant directories with abundance files

### merge
- Purpose: aggregate per-sample quantifications into a single matrix
- Common params: `out` (Path to merged table, default `work_dir/merged_abundance.tsv`)
- Artifacts: merged counts/TPMs table suitable for downstream normalization/DE

### cstmm
- Purpose: normalization/scaling of merged counts
- Common params: `out-dir` (default `work_dir/cstmm`)
- Artifacts: normalized matrices and summaries

### curate
- Purpose: optional data curation/cleanup steps (samples/genes), per amalgkit
- Artifacts: curated datasets and reports under `work_dir/curate`

### csca
- Purpose: clustering/sample assessment and QC visualization
- Common params: `out-dir` (default `work_dir/csca`)
- Artifacts: clustering results, plots, metrics

### sanity
- Purpose: final checks to ensure a valid, complete run
- Artifacts: status files and logs

## Downstream meta-analysis

Use the merged (and optionally normalized) matrices for:
- Within-study DE: DESeq2/edgeR
- Cross-study meta-analysis: combine per-gene effects/p-values (e.g., metaRNASeq). Batch-effect handling can leverage ComBat/limma on matrices prior to DE or within model frameworks.

All files should reside under `output/` by default per repository policy. Use `work_dir` to relocate artifacts.
