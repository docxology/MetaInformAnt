# Agent Directives: docs/rna/amalgkit

## Role
Documentation for amalgkit RNA-seq workflow integration.

## Directory Structure
- `steps/` - Step-by-step documentation for the nine-stage per-species chain

## Key Files
- `amalgkit.md` - Complete amalgkit integration guide
- `FUNCTIONS.md` - API function reference
- `commands.md` - CLI command reference
- `genome_preparation.md` - Genome setup for quantification
- `monitoring.md` - SQLite, log, and checkpoint status inspection
- `R_INSTALLATION.md` - optional downstream R environment setup

## Workflow steps (in `steps/`)

1. `metadata` - Fetch sample metadata
2. `select` - Apply the reviewed selection rules
3. `getfastq` - Acquire and validate reads
4. `integrate` - Attach validated local read paths
5. `quant` - Quantify with Kallisto
6. `merge` - Merge current sample quantifications
7. `wsfilter` - Within-species filtering
8. `finalize` - Produce analysis-ready tables
9. `sanity` - Validate output and provenance

`cstmm` and `csfilter` are opt-in cross-species commands and are documented
separately from the per-species chain.
