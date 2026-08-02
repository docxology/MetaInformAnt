# Agent Directives: docs/rna

## Role

Documentation for the RNA transcriptomics module.

## Module Scope

ENA-first amalgkit streaming pipeline, Kallisto quantification, cross-species
TMM normalization, and preemption-aware orchestration. Sample totals are
data-root dependent; historical production counts must not be presented as
current inventory without a generated evidence report.

## Key Source Files

| Path | Description |
|------|-------------|
| `src/metainformant/rna/amalgkit/` | Amalgkit CLI integration |
| `src/metainformant/rna/engine/` | Workflow orchestration |
| `src/metainformant/rna/retrieval/` | ENA/SRA data retrieval |
| `src/metainformant/rna/analysis/` | Expression matrix analysis |
