# Personal AI Infrastructure (PAI) - docs/rna/amalgkit/steps

## Context & Intent

- **Path**: `docs/rna/amalgkit/steps/` (relative to repo root)
- **Purpose**: Per-command reference pages for the Amalgkit 0.16.60 workflow:
  inputs, outputs, options, and common issues per numbered stage.
- **Domain**: `rna/amalgkit`. CLI source: pinned Amalgkit 0.16.60; Python
  wrappers in `src/metainformant/rna/`.

## Virtual Hierarchy

- **Type**: Documentation
- **Parent**: `docs/rna/amalgkit/`
- **Children**: none. 18 files: 14 step/topic pages: 01_metadata.md, 02_dataset.md, 03_select.md, 04_getfastq.md,
  05_integrate.md, 06_quant.md, 06_quant_advanced.md, 06_quant_troubleshooting.md, 07_merge.md, 08_cstmm.md,
  09_wsfilter.md, 10_csfilter.md, 10_finalize.md, 11_sanity.md — plus AGENTS/PAI/README/SPEC.

## Maintenance Notes

From the sibling `AGENTS.md` / `README.md` / `SPEC.md`:

- Read in order; each file documents inputs, outputs, options, issues
  (`AGENTS.md`); the CLI is authoritative — verify `amalgkit <cmd> --help`
  before documenting options (`README.md`).
- `02_dataset` is workspace preparation, not a biological stage; the chain is
  metadata -> select -> getfastq -> integrate -> quant -> merge -> wsfilter ->
  finalize -> sanity (`README.md`).
- `08_cstmm` and `10_csfilter` are optional cross-species branches; a
  non-empty directory is not evidence of completion; tables tie to
  metadata + configuration hash (`README.md`).

## AI Workflows

- **Diagnose a failing stage**: open the matching numbered page (e.g.
  `06_quant_troubleshooting.md`) before touching `src/metainformant/rna/`.
- **Document a new option**: confirm against 0.16.60 help output, then update
  the numbered page and the table in `README.md`.
- **Validation**: repo docs checks (`scripts/verify_documentation_code.py`,
  internal-link validation); no dir-local test suite.
