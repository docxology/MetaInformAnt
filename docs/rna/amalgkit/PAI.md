# Personal AI Infrastructure (PAI) - docs/rna/amalgkit

## Context & Intent

- **Path**: `docs/rna/amalgkit/` (relative to repo root)
- **Purpose**: Documentation for the Amalgkit 0.16.60 cross-species RNA-seq
  integration: command chain, genome prep, ortholog generation, monitoring.
- **Domain**: `rna` — implementation in `src/metainformant/rna/` + `projects/hymenoptera_amalgkit/`.

## Virtual Hierarchy

- **Type**: Documentation
- **Parent**: `docs/rna/`
- **Children**: `steps/` (per-command reference, 01-11 numbered pages).
- **This directory**: 17 topic files besides `steps/`: AGENTS.md, FUNCTIONS.md,
  PATH_RESOLUTION.md, README.md, R_INSTALLATION.md, SPEC.md, TROUBLESHOOTING.md,
  amalgkit.md, commands.md, cross_species_pipeline.md, genome_preparation.md,
  genome_setup_guide.md, guide.md, monitoring.md, ortholog_generation.md,
  testing_coverage.md, tissue_patching.md.

## Maintenance Notes

From the sibling `AGENTS.md` / `README.md`:

- CLI authoritative: "verify help output before adding a command-specific
  option"; "`cstmm` and `csfilter` are opt-in cross-species commands"
  (`AGENTS.md`); real-data rules (`steps/README.md`): a non-empty directory
  is not evidence of completion; tables tie to metadata + config hash.
- Evidence chain: producer records SQLite sample state, hashes, and resumable
  downloads; the downstream runner verifies receipts before reusing matrices.

## AI Workflows

- **Understand the chain**: `README.md` -> `commands.md` -> `steps/` pages.
- **Run**: `bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh`
  (full) or `scripts/rna/process_species.py --species apis_mellifera ...`
  (diagnostic); never run while the live producer holds locks.
- **Troubleshoot**: `TROUBLESHOOTING.md` and `PATH_RESOLUTION.md` first.
- **Validation**: `uv run python scripts/rna/validate_configs.py` — no
  dir-local test suite.
