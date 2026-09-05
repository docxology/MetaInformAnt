# Personal AI Infrastructure (PAI) - docs/singlecell/config

## Context & Intent

- **Path**: `docs/singlecell/config/` (relative to repo root)
- **Purpose**: Documentation directory reserved for single-cell RNA-seq
  configuration guidance (filtering, normalization, clustering resolution,
  dimensionality reduction, integration), per the sibling `AGENTS.md`.
- **Domain**: `singlecell`. Real implementation:
  `src/metainformant/singlecell/`; parent topic docs in `docs/singlecell/`.

## Virtual Hierarchy

- **Type**: Documentation
- **Parent**: `docs/singlecell/` (topics: preprocessing, clustering,
  dimensionality, integration, celltyping, differential, trajectory,
  velocity, visualization — .md).
- **Children**: none. Exactly four files, all markdown: `AGENTS.md`, `PAI.md`,
  `README.md`, `SPEC.md`. **No example config files exist here yet** — the
  README states "No runnable files are stored directly in this directory yet."

## Maintenance Notes

From the sibling `AGENTS.md`:

- "Reference these configs when setting up single-cell analysis workflows"
  describes intended scope; until example files exist, keep claims limited
  to what exists.
- Sibling `README.md` rule: keep this reference synchronized with the
  corresponding source module and domain guide (`src/metainformant/singlecell/`,
  `docs/singlecell/`).

## AI Workflows

- **Configuration work today**: source actual knobs from `src/metainformant/singlecell/` code.
- **Adding examples**: add real YAML/config files copied from runnable
  analyses (real-implementation policy: no placeholder configs), then update
  `README.md` and `SPEC.md` in the same change.
- **Validation**: repo docs checks (`scripts/verify_documentation_code.py`,
  internal-link validation); no dir-local test suite exists.
