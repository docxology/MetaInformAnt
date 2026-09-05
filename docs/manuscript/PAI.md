# Personal AI Infrastructure (PAI) - docs/manuscript

## Context & Intent

- **Path**: `docs/manuscript/` (relative to repo root)
- **Purpose**: Docs-only directory recording that the parent package has
  **no publication-track manuscript** and what would trigger creating one;
  the single substantive artifact is `MANUSCRIPT_STATUS.md` (2026-09-05).
- **Domain**: `docs` — there is no corresponding `src/metainformant/` module
  and no `manuscript/` tree at the repo top level.

## Virtual Hierarchy

- **Type**: Documentation
- **Parent**: `docs/` (doc-tree entry point: `docs/index.md`)
- **Children**: none. Flat directory of exactly five files: `AGENTS.md`,
  `MANUSCRIPT_STATUS.md`, `PAI.md`, `README.md`, `SPEC.md`.

## Maintenance Notes

Quoted from the sibling `AGENTS.md` (verified 2026-08-29):

- "Docs-only directory — no corresponding `src/metainformant/` module."
- "Keep prose aligned with the current checkout, not aspiration."
- "Cross-link to existing local files only."

- Per `MANUSCRIPT_STATUS.md`, a real manuscript tree (config.yaml, section
  files 00-99, references.bib) belongs at the **repo top level**, not here;
  per-study manuscripts belong to nested `projects/` subprojects (own
  branch/PR).

## AI Workflows

- **Before any manuscript request**: read `MANUSCRIPT_STATUS.md`; do not
  invent variables, figures, or citation ledgers —
  `docs/reviews/RESEARCH_SOFTWARE_REVIEW_2026-08-13.md` withholds those.
- **Update**: revise `MANUSCRIPT_STATUS.md` only when repo structure changes
  (a subproject manuscript appears, or a methods paper starts).
- **Validation**: run repo docs checks (`scripts/verify_documentation_code.py`,
  internal-link validation) after edits; there is no dir-local test.
