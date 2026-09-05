# Specification: docs/manuscript

## Scope

- Documentation for the "manuscript" area of METAINFORMANT: exactly five
  markdown files — `AGENTS.md`, `MANUSCRIPT_STATUS.md`, `PAI.md`,
  `README.md`, `SPEC.md` (verified 2026-09-05; no hidden subdirectories).
- Documents the verified absence of a parent-package manuscript: the repo is
  an analysis platform whose deliverables are code, pipelines, configs, and
  generated analyses, so no publication target applies today.

## Architecture

- Docs-only, flat, no code and no build step; no generated artifacts.
- Governance: sibling `AGENTS.md` rules plus the repository-root `AGENTS.md`.
- Boundary: per-study manuscripts live in nested subprojects
  (`projects/apis_gwas`, `projects/hymenoptera_amalgkit`,
  `projects/drosophila_scrna_2026`) and stay inside those repository
  boundaries, governed by their own PRs.

## Data Structures

- No Python modules, no config files, no `references.bib`; prose markdown only.
- Key artifact `MANUSCRIPT_STATUS.md` is a status ledger with sections:
  repo type, evidence checked, why no manuscript applies, and trigger
  conditions — a methods paper for the platform itself, or promotion of one
  subproject's findings — at which point a full `manuscript/` tree
  (config.yaml, section files 00-99, references.bib) is created at repo top
  level following the docxology/template standard.

## API Definition

- No API. Consumers read `MANUSCRIPT_STATUS.md` directly.
- Cross-references (existing local files only): `docs/index.md`,
  `docs/CAPABILITY_MATRIX.md`,
  `docs/reviews/RESEARCH_SOFTWARE_REVIEW_2026-08-13.md`.
- Checks: covered by repo documentation validation
  (`scripts/verify_documentation_code.py`, internal-link validation);
  no dir-local test exists.
