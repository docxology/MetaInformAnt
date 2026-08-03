# RNA module specification

## Scope

Version-pinned RNA-seq acquisition, quantification, evidence validation, and
downstream analysis for the current Amalgkit 0.16.33 contract.

## Current packages

| Package | Responsibility |
|---|---|
| `amalgkit/` | Version checks, command registry, and CLI wrappers |
| `analysis/` | Expression, QC, cross-species, and validation calculations |
| `core/` | Configuration, cleanup, dependencies, and environment checks |
| `engine/` | Streaming production, workflow planning, SQLite progress, and provenance |
| `retrieval/` | ENA/SRA retrieval helpers |
| `splicing/` | Alternative-splicing analysis |

## Current stages

The production contract is:

`metadata → select → getfastq → integrate → quant → merge → wsfilter → finalize → sanity`

The first five stages are resumable sample production. The final four stages
are project checkpoints and require the producer lock to be released.

## Public interfaces

- `engine.streaming_orchestrator.StreamingPipelineOrchestrator` — bounded
  multi-species producer with ENA-first acquisition and resumable quantification
- `engine.progress_db.ProgressDB` — SQLite state machine for sample progress
- `engine.workflow` — `AmalgkitWorkflowConfig`, `load_workflow_config`,
  `plan_workflow`, and `execute_workflow`
- `engine.provenance` — current hash-bound receipt writers and validators
- `amalgkit` — current Amalgkit version and command registry
- `analysis` — public expression, QC, and cross-species namespaces

## Invariants

- Species configuration is selected from the project config directory.
- The data root is explicit and external to Git.
- Repeated discovery and dry-run operations are read-only and deterministic.
- Every resumable sample state is stored in SQLite.
- Quantification reuse requires the current provenance sidecar and exact file
  digest.
- Downstream matrices are admitted only after current metadata, selection,
  quantification, and stage receipts validate.
- Biological claims are made only from regenerated, source-bound evidence.
