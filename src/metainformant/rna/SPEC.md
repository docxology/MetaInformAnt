# Specification: rna

## 🎯 Scope

RNA transcriptomic analysis and workflow orchestration module for METAINFORMANT.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## 📦 Sub-packages

| Sub-package | Purpose |
|---|---|
| `amalgkit/` | Amalgkit CLI wrapper, genome prep, metadata filtering, tissue normalization |
| `analysis/` | Expression analysis, QC metrics/filtering, cross-species, validation |
| `core/` | Configs, cleanup, dependency checks, environment setup |
| `deconvolution/` | Cell-type proportion estimation from bulk RNA-seq |
| `engine/` | Workflow planning/execution, orchestration, monitoring, discovery |
| `retrieval/` | ENA/SRA data download |
| `splicing/` | Alternative splicing detection, isoform quantification |

## 🔌 API Definition

### Key Exports

- `engine.workflow` — Re-export hub: `load_workflow_config`, `execute_workflow`, `plan_workflow`
- `engine.streaming_orchestrator.StreamingPipelineOrchestrator` — Multi-species ENA-first pipeline
- `engine.orchestrator.StreamingPipeline` — Per-sample download→quant→cleanup pipeline
- `amalgkit.amalgkit` — CLI wrapper: `run_amalgkit`, `build_amalgkit_command`, step functions
- `amalgkit.AmalgkitParams` — Typed parameter container
- `analysis.validation` — Sample/pipeline validation
- `retrieval.ena_downloader.ENADownloader` — Direct ENA FASTQ downloader
