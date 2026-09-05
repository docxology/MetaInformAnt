# Specification: docs/reviews

## Scope

- Documentation for the "reviews" area of METAINFORMANT: dated review
  ledgers over the parent package. Exactly five markdown files — `AGENTS.md`,
  `PAI.md`, `README.md`, `RESEARCH_SOFTWARE_REVIEW_2026-08-13.md`, `SPEC.md`
  (verified 2026-09-05).
- `RESEARCH_SOFTWARE_REVIEW_2026-08-13.md` records the verified outcome
  "INCOMPLETE", nine review lanes (parent architecture passed; MCP
  scaffold/not applicable; parent manuscript pipeline not applicable;
  biological inference blocked, etc.), implemented security improvements
  (hardened NCBI ZIP extraction replacing `ZipFile.extractall()`), and a
  validation ledger (e.g. full pytest 7,708 passed / 27 skipped;
  `scripts/core_docs_cross_check.py --strict` 235 symbols, 0 issues;
  3,269 internal links, 0 broken; mypy ratchet 169 errors within budget 171).

## Architecture

- Docs-only, flat, append-only ledger style; no code and no build step.
- Governance: sibling `AGENTS.md` rules plus the repository-root `AGENTS.md`.
- Boundary: nested projects (`projects/apis_gwas` BeeWAS, Hymenoptera
  campaign) are reviewed but governed inside their own repositories.

## Data Structures

- No Python modules or configs; prose markdown only.
- Review files use `##` sections: Verified outcome; Scope and applicability
  (lane/status/evidence table); Improvements implemented; Current source and
  generated surfaces; Validation ledger (command/status/result table);
  Scientific boundaries; Remaining gates and next actions.

## API Definition

- No API. Consumers read the dated review files directly.
- Cross-references (existing local files only): `docs/index.md`,
  `docs/CAPABILITY_MATRIX.md`, `docs/MIGRATION_0.4.md`,
  `docs/tasks/VALIDATION_REPORT.md` (historical), `docs/manuscript/MANUSCRIPT_STATUS.md`.
- Checks: repo documentation validation (`scripts/verify_documentation_code.py`,
  internal-link validation); no dir-local test exists.
