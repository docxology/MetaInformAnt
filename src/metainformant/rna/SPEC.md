# Specification: rna

## Scope

RNA transcriptomic analysis and workflow orchestration module for METAINFORMANT.

## Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## Sub-packages

| Sub-package | Purpose |
|---|---|
| `amalgkit/` | Amalgkit CLI wrapper, genome prep, metadata filtering, tissue normalization |
| `analysis/` | Expression analysis, QC metrics/filtering, cross-species, validation |
| `core/` | Configs, cleanup, dependency checks, environment setup |
| `deconvolution/` | Cell-type proportion estimation from bulk RNA-seq |
| `engine/` | Workflow planning/execution, orchestration, monitoring, discovery |
| `retrieval/` | ENA/SRA data download |
| `splicing/` | Alternative splicing detection, isoform quantification |

## API Definition

### Key Exports

- `engine.workflow` — Re-export hub: `load_workflow_config`, `execute_workflow`, `plan_workflow`
- `engine.streaming_orchestrator.StreamingPipelineOrchestrator` — Multi-species ENA-first pipeline
- `engine.orchestrator.StreamingPipeline` — Per-sample download→quant→cleanup pipeline
- `amalgkit.amalgkit` — CLI wrapper: `run_amalgkit`, `build_amalgkit_command`, step functions
- `amalgkit.AmalgkitParams` — Typed parameter container
- `analysis.validation` — Sample/pipeline validation
- `analysis.expression_core` — Count normalization (`normalize_counts`), size factors, filtering, highly variable genes
- `analysis.expression_analysis` — Differential expression, p-value adjustment, PCA, distances, volcano/MA table preparation
- `analysis.expression` — Backward-compatible shim re-exporting `expression_core` and `expression_analysis`
- `analysis.qc_metrics` — Numeric QC metrics for samples/genes, outliers, library complexity, saturation, correlation matrices
- `analysis.qc_filtering` — Batch effect detection, GC/length bias detection, and QC report generation
- `analysis.qc` — Backward-compatible shim re-exporting QC metrics/filtering helpers
- `analysis.cross_species` — Ortholog mapping, conservation scoring, divergence matrices, phylogenetic profiles, cross-species PCA
- `analysis.protein_integration` — RNA/protein translation efficiency (`ratio`, `correlation`), linear protein-abundance prediction, ribosome profiling integration
- `retrieval.ena_downloader.ENADownloader` — Direct ENA FASTQ downloader

### Analysis Validation Contracts

- QC count matrices are numeric, finite, and non-negative; empty matrices raise `ValueError`.
- Batch labels must include all expression samples; extra labels are ignored.
- GC content for matched genes must be finite and in `[0, 1]`; matched gene lengths must be positive.
- RNA/protein integration rejects unsupported methods explicitly and filters NaN measurements deterministically.
