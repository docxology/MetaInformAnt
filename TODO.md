# MetaInformAnt TODO

> **Last Updated**: 2026-08-09
>
> This TODO contains only active unfinished work. Prior BeeWAS GWAS,
> manuscript, dashboard, validation, and release-bundle work is recorded in
> status notes and generated release artifacts, not in this active list.

The Hymenoptera Amalgkit project has a separate, data-root-scoped backlog at
[`projects/hymenoptera_amalgkit/TODO.md`](projects/hymenoptera_amalgkit/TODO.md).
Keep that work there rather than mixing a live external-volume transcriptome
campaign with the BeeWAS release backlog below. Completed MetaInformAnt RNA
and Hymenoptera implementation work is intentionally absent from both active
lists; its evidence is maintained in the project status and readiness records.

---

## Current Scope

BeeWAS four-population statistical GWAS reporting is operational for the
combined cohort and population scopes `C`, `I`, `M`, and `R`. Manuscript
rendering, visual QA, artifact indexing, guarded dashboard generation, GO
curation-atlas outputs, statistical deep-dive summaries, and release packaging
are generated and validated. The active work below is evidence-gated: it
requires human-approved metadata, human-reviewed GO mapping decisions, or
independent replication evidence before the repository can mark the item
complete.

Use the package-module sweep to refresh status after any evidence is supplied:

```bash
cd /Users/mini/Documents/GitHub/MetaInformAntCode/projects/apis_gwas
uv run python -m beewas.gwas.review_release build
```

---

## Package-wide review backlog (2026-08-09)

This section records the package-level improvement scope from the stewardship
review. It is a review backlog, not a claim that these items are fixed.
Preserve existing uncommitted parent work, nested repository boundaries, and
the live RNA campaign while addressing it. Historical audit reports, negative
controls, failure evidence, and provenance artifacts remain part of the record.

### P1 — Correctness, data safety, and scientific readiness

- [x] **Make the NCBI contact contract explicit.** Explicit contact values take
  precedence over `NCBI_EMAIL`; otherwise callers must opt into anonymous
  access. Anonymous mode warns, leaves `Entrez.email` unset, and does not
  persist invented identity data. Offline tests cover each mode.
  Scope: `src/metainformant/rna/engine/discovery.py`,
  `src/metainformant/gwas/data/sra_download.py`,
  `src/metainformant/dna/external/genomes.py`, and shared configuration/tests.

- [x] **Make SRA configuration campaign-local and non-destructive.** Remove
  production dependence on `vdb-config -s` against user-global NCBI
  configuration and stop cleanup routines from implicitly scanning or moving
  files under `Path.home()/ncbi` or `/tmp/ncbi`. Native Amalgkit, direct
  `fasterq-dump`, validation subprocesses, and cleanup must inherit one
  explicit campaign-local environment and path contract.
  Scope: `src/metainformant/rna/amalgkit/sra_environment.py`,
  `src/metainformant/rna/engine/workflow_steps.py`,
  `src/metainformant/rna/engine/workflow_cleanup.py`, and subprocess tests.

- [x] **Make readiness reports provenance- and manifest-driven.** Replace file
  existence and `pgrep` completion heuristics in
  `scripts/rna/check_pipeline_status.py` with lock ownership, durable
  receipts, current manifests, and producer heartbeats. Report executable
  readiness, cohort readiness, descriptive analysis, and biological inference
  as separate states; never promote a partial cohort or descriptive matrix to
  a biological result.
  Exit evidence: fixture reports for empty, active, partial, stale, current,
  and finalized states, each independently auditable.

- [x] **Separate quantification reuse from cleanup authorization.** A
  version-drift or otherwise compatible output may be reusable for inspection,
  but raw deletion must require exact current quantification provenance plus a
  current manifest and receipt. Preserve failed, partial, stale, legacy, and
  version-drift states as visible/resumable records.
  Scope: `src/metainformant/rna/engine/provenance.py`,
  `src/metainformant/rna/engine/progress_db.py`,
  `scripts/rna/reclaim_quantified_raw.py`, and orchestrator reset/retry policy.

- [x] **Make provenance identity independent of the caller's working
  directory.** Canonicalize configuration and output paths at receipt
  creation/classification boundaries and test classification from a different
  process CWD. Keep hash-bound output reuse strict and observable.

- [x] **Require explicit cross-species sample/condition alignment.** Remove or
  quarantine positional column fallback for conservation/divergence comparisons
  unless the caller supplies a validated alignment map and metadata evidence.
  Define duplicate-label, missingness, replicate, and condition semantics
  before any comparison is called conserved or divergent.
  Scope: `src/metainformant/rna/analysis/cross_species.py`, visualization
  consumers, RNA docs, and synthetic/fixture tests.

- [ ] **Make cross-species statistics honest and reproducible.** Specify
  estimands and replicate units, apply an appropriate multiple-testing
  procedure, distinguish descriptive permutation scores from inferential
  p-values, and record random seeds/parameters in provenance. Add tree and
  orthology invariants for phylogenetic profiles.

### P1 — Package quality and release gates

- [x] **Repair and then enforce formatting/lint quality across maintained
  code.** The current changed-tree review found Black issues in six changed
  files, isort issues in three, and flake8 issues including E402/E501, an
  unused import, and a nested-definition spacing error. Include maintained
  `scripts/` in the blocking gate while keeping compatibility/generated
  exceptions narrow and documented.
  Exit evidence: check-only Black/isort/flake8 passes for `src/`, `tests/`, and
  maintained scripts; CI fails on new findings rather than masking them.

- [x] **Make typing claims and gates truthful.** `SPEC.md` and contributor
  guidance imply strict typing, while the current RNA mypy review reports 173
  errors across 42 files and the configured gate is not strict. Choose a
  phased baseline: type core/RNA contracts first, add required third-party
  stubs or narrow facades, and enforce only the documented scope at each phase.
  Exit evidence: a checked-in typing policy, a ratcheting error baseline or
  clean scoped gate, and tests for the canonical core, DNA, RNA, and NCBI
  contracts. Removed compatibility paths are documented only in
  `docs/MIGRATION_0.4.md`.

- [x] **Remove non-blocking quality semantics where release correctness
  depends on them.** Reassess `continue-on-error` for extended/network/all/FAT/
  example jobs, optional documentation builds, masked report generation, and
  package validation that only warns when required tooling is absent. Keep
  genuinely optional jobs explicitly labeled informational.
  Scope: `.github/workflows/test.yml`, `.github/workflows/release.yml`,
  `scripts/package/test.sh`, `scripts/package/uv_quality.sh`, and
  `scripts/package/build.sh`.

- [x] **Add disposable installation smoke tests to release validation.** Build
  the wheel, install it into clean Python 3.11 and 3.12 environments, test
  public and legacy imports, verify the CLI/API surface, and run documented
  lightweight checks. Do not recreate or mutate the live producer environment.

### P2 — Documentation, generated surfaces, and capability truth

- [x] **Repair command/documentation drift.** Replace references to the
  retired test entry point with the maintained package test entry point, then
  audit active README, contributor, CI, scripts, and RNA docs for
  retired commands and paths. Keep historical reports explicitly historical.
  Scope: `CONTRIBUTING.md`, `docs/CONTRIBUTING.md`, `tests/README_tests.md`,
  active README/docs, and release instructions.

- [x] **Publish a capability matrix for the package.** Classify each public
  module/API as stable, experimental, scaffold/placeholder, deprecated
  compatibility facade, or retired. Reconcile active README claims with
  validation reports documenting placeholders or missing implementations,
  especially tasks, life-events, networks, ML, math, and variant-download
  surfaces.
  Exit evidence: a source-linked matrix, promotion criteria, and tests/docs
  that do not imply scientific or operational completeness for scaffolds.

- [x] **Keep generated API documentation and Cursor skills reproducible.** Make
  generation/check commands part of the required documentation gate, verify
  internal links, and ensure generated outputs record source/version
  provenance without overwriting historical evidence. Current review evidence:
  core cross-check, Cursor skill check, RNA docs validation, and internal link
  validation pass.

- [x] **Define package version and release authority once.** Reconcile the
  package version in `pyproject.toml`, `src/metainformant/__init__.py`,
  generated metadata, changelog/release notes, and hosted release workflow.
  Document which generated artifacts must be refreshed and which approvals are
  required before publication.

### P2 — Tests, contracts, and security

- [x] **Expand contract tests around the current RNA boundary.** The current
  targeted baseline covers
  fixtures for lock mutual exclusion, interrupted recovery, hash-bound reuse,
  provenance migration, raw-cleanup refusal, subprocess environments, alias
  resolution, partial cohorts, and downstream finalization.

- [x] **Add package-wide public-surface and optional-dependency tests.** Exercise
  maintained module imports, legacy aliases, a clean wheel install, no-network
  defaults, real-implementation policy, and documented optional dependency
  behavior. Structural completeness and export audits are useful guards but do
  not establish scientific correctness; pair them with behavior-level
  contracts.

- [ ] **Run blocking security and dependency checks on the release surface.**
  Review Bandit/Safety (or selected maintained equivalents), lockfile
  provenance, subprocess argument boundaries, path traversal, temporary-file
  handling, external-tool logs, and contact-data persistence. Record justified
  compatibility exceptions rather than masking findings.

### P3 — Operations and downstream analysis

- [ ] **Characterize campaign-scale performance without touching the live
  producer.** Benchmark discovery, bounded acquisition/quantification queues,
  retry/backoff, hash/provenance I/O, storage pressure, and resource limits in
  disposable fixtures. Define budgets and telemetry for long-running
  campaigns; never start a competing producer during validation.

- [ ] **Make evidence bundles first-class release artifacts.** Include immutable
  input/config hashes, software/tool versions, lock/heartbeat history, terminal
  and unresolved task counts, receipts, manifests, matrices, error taxonomy,
  and command transcripts. State acquisition, recovery, quantification,
  descriptive analysis, and biological inference separately.

- [ ] **Specify downstream analytical contracts before promotion.** Document
  merge/sample/ortholog identifiers, normalization, differential-expression or
  equivalent methods, batch handling, replicate requirements, missing-data
  policy, statistical corrections, and the evidence required to move from
  descriptive analysis to biological inference.

### Boundary and execution rules

- Do not edit nested repository contents from the parent. Reconcile nested
  gitlinks only after their own branches/PRs contain validated, remotely
  reachable commits.
- Do not touch `/Volumes/blue/data/amalgkit` except for explicitly read-only
  monitoring.
- Do not reset, discard, recreate the active environment, or claim
  campaign/scientific completion from a moving operational snapshot.
- Before each release claim, recapture parent/nested commits, reachability,
  dirty paths, PR state, campaign status, current manifest, provenance, and
  reproducible evidence-bundle status.

---

## Minor Updates

- [ ] **Finalize BeeWAS manuscript metadata** - Complete
  `projects/apis_gwas/doc/manuscript/renderable/final_metadata_intake.tsv`
  with final title, authors, affiliations, corresponding author, funding,
  contributions, conflicts, data availability, code availability, reviewer
  name, review date, and final status values. Current status is tracked in
  `projects/apis_gwas/output/pdf/beewas_gwas/manuscript_metadata_status.tsv`.

## Medium Improvements

- [ ] **Review BeeWAS GO mapping evidence** - Review
  `projects/apis_gwas/doc/operations/reviews/gwas_reviewed_go_mapping_intake_template.tsv`,
  set justified rows to `review_status=reviewed`, retain provenance fields, and
  ingest the reviewed TSV with `uv run python -m beewas.gwas.reviewed_go ingest
  path/to/reviewed_go.tsv`. Current intake status is tracked in
  `projects/apis_gwas/results/gwas_79cram_clean_reviewed_fullgemma_guarded_20260624_summary/tables/gwas_reviewed_go_mapping_intake_status.tsv`.

## Major Improvements

- [ ] **Supply independent replication evidence** - Complete
  `projects/apis_gwas/doc/operations/reviews/gwas_replication_review_scaffold.tsv`
  with external cohort details, comparable trait definitions, model
  compatibility, evidence URIs, reviewer decisions, reviewer names, and review
  dates before any replication claim is made.

---

## Operating Rules

- Use package-module commands for BeeWAS work, for example
  `uv run python -m beewas.gwas.guarded ...`.
- Keep runtime outputs under `output/`, `results/`, or documented data/output
  directories.
- For the live Hymenoptera cloud campaign, reattach an interrupted controller
  with `--skip-bundle-upload` and the original state/source/input IDs; never
  rebuild or overwrite an immutable source object under an existing run ID.
- Do not convert curation evidence into biological, causal, gene-function,
  replication, breeding, or final heritability claims without a reviewed
  decision artifact.
- Keep this TODO forward-looking; move resolved evidence and dated snapshots to
  the appropriate status or release document.
