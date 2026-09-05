# Personal AI Infrastructure (PAI) - docs/reviews

## Context & Intent

- **Path**: `docs/reviews/` (relative to repo root)
- **Purpose**: Documentation-only directory holding dated review ledgers.
  Current sole artifact: `RESEARCH_SOFTWARE_REVIEW_2026-08-13.md`, whose
  verified outcome is "INCOMPLETE — engineering hardening advanced; release
  and scientific promotion remain gated."
- **Domain**: `docs` — no corresponding `src/metainformant/` module.

## Virtual Hierarchy

- **Type**: Documentation
- **Parent**: `docs/` (doc-tree entry point: `docs/index.md`)
- **Children**: none. Flat directory of exactly five files: `AGENTS.md`,
  `PAI.md`, `README.md`, `RESEARCH_SOFTWARE_REVIEW_2026-08-13.md`, `SPEC.md`.

## Maintenance Notes

Quoted from the sibling `AGENTS.md` (verified 2026-08-29):

- "Docs-only directory — no corresponding `src/metainformant/` module."
- "Keep prose aligned with the current checkout, not aspiration."
- "Cross-link to existing local files only."

- The review file states: "This record is a review ledger, not a scientific
  result." Never promote its lane statuses into claims; it withholds
  manuscript, statistical-validity, and biological-inference statements.

## AI Workflows

- **Before citing release readiness**: read the review's "Scope and
  applicability" and "Validation ledger" tables; cite the dated file, never
  paraphrase statuses as current.
- **New review**: add a date-stamped file
  (`RESEARCH_SOFTWARE_REVIEW_<YYYY-MM-DD>.md`); never edit historical
  ledgers — they stay labeled and separate from current evidence.
- **Validation**: run repo docs checks (`scripts/verify_documentation_code.py`,
  internal-link validation) after edits; no dir-local test exists.
