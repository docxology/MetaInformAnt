# Research-software and scholarship review — 2026-08-13

## Verified outcome

**INCOMPLETE — engineering hardening advanced; release and scientific promotion
remain gated.** The parent package has a current source/test/documentation
contract on branch `agent/publish-metainformant-full-20260802`; the preceding
parent head is published to draft PR #192 and this review candidate is being
validated before its update. The live RNA/Amalgkit producer was observed active
during this review and was not interrupted. No cohort-complete,
biological-inference, manuscript-publication, or deposit claim is made.

This record is a review ledger, not a scientific result. It intentionally keeps
historical audit reports and nested-repository evidence separate from current
source-bound release evidence.

## Scope and applicability

| Review lane | Status | Evidence and boundary |
| --- | --- | --- |
| Parent package architecture and contracts | passed | `src/metainformant/`, `tests/`, package verification, and generated API checks |
| External-tool and path security | passed for reviewed boundaries | Campaign-local SRA environments, accession validation, and confined NCBI ZIP extraction |
| RNA operational readiness | available but not complete | Current locks, receipts, manifests, and producer heartbeat govern promotion; active acquisition remains a hard gate |
| GWAS download and SRA contracts | passed for local contracts | Offline-safe subprocess environment and accession/path validation are tested; live NCBI availability remains external |
| MCP | scaffold / not applicable as a server release | Only the standalone Amalgkit monitor is implemented; no server, registry, or transport is documented as available |
| Parent manuscript/render pipeline | not applicable | The parent package has no canonical manuscript variable/render/figure registry pipeline |
| BeeWAS manuscript and scholarship | nested / separate owner | Review source exists under `projects/apis_gwas`; its own branch, tests, provenance, and draft PR govern changes |
| Hymenoptera campaign content | nested / preserved | The nested checkout has pre-existing dirty work and was not edited from the parent |
| Biological inference | blocked | Producer is active and current cohort/manifests/evidence are not a completed scientific result |

## Improvements implemented in the current review pass

- Replaced unrestricted `ZipFile.extractall()` in both NCBI Datasets genome and
  annotation paths with explicit extraction that rejects absolute paths,
  traversal components, resolved destination escapes, and symlink members.
- Added adversarial regression tests for traversal rejection, symlink rejection,
  and valid nested extraction.
- Clarified the pre-commit contract so repository-wide mypy remains a ratcheted
  verification gate and documentation style is explicitly informational while
  blocking link/API checks remain authoritative.
- Preserved the shared campaign-local SRA environment contract across RNA and
  GWAS subprocesses; no user-global `vdb-config` or home-directory cleanup was
  introduced.

## Current source and generated surfaces

- API source inventory: 235 public symbols across 25 modules; strict core-doc
  cross-check reported zero issues.
- Cursor skills: 167 generated skills match their `AGENTS.md` locations.
- Internal documentation links: 3,269 checked and zero broken in the latest
  completed link-validation pass.
- Capability and migration contracts: `docs/CAPABILITY_MATRIX.md` and
  `docs/MIGRATION_0.4.md` are the source-linked 0.4.0 release boundary.
- Historical audit material remains labeled and is not treated as current
  evidence; for example, `docs/tasks/VALIDATION_REPORT.md` explicitly identifies
  its 2026-04-29 snapshot as historical.

The parent has no canonical manuscript variable registry, source-to-render
hydration pipeline, figure registry, PDF renderer, or bibliography ledger. The
review therefore does not fabricate those artifacts. The nested BeeWAS project
does contain such manuscript scaffolding, but its contents remain within the
nested repository boundary and are governed by its independent PR.

## Validation ledger

Statuses use the review brief’s vocabulary. Results below are from the current
candidate before this ledger commit, except where a command is explicitly
marked pending.

| Command or check | Status | Result / limitation |
| --- | --- | --- |
| `env -u VIRTUAL_ENV uv run pytest -q` | passed | 7,708 passed, 27 skipped |
| Focused GWAS SRA/security tests | passed | Existing SRA tests plus new ZIP adversarial tests pass |
| Focused RNA environment/cleanup tests | passed | 51 passed in the latest focused run |
| Black, isort, flake8 on changed files | passed | Blocking checks pass after formatting |
| Bandit `-ll` on changed files | passed | No medium/high findings |
| `scripts/core_docs_cross_check.py --strict` | passed | 235 symbols, 0 issues |
| `scripts/package/generate_cursor_skills.py --check` | passed | 167 skills match source locations |
| RNA configuration validation | passed | 27 species configurations; Amalgkit 0.16.60 |
| Mypy ratchet | passed with budget | 169 errors in 41 files, within the documented budget of 171; package is not type-clean |
| Internal link validation | passed | 3,269 links checked, 0 broken |
| Disposable Python 3.11/3.12 wheel installs | not run in this pass | Hosted matrix and package smoke coverage remain release evidence to refresh |
| Hosted PR checks after latest push | pending | New checks were in progress when this record was prepared |
| Live campaign completion | blocked | Active producer and unresolved acquisition state; no interruption performed |

## Scientific, manuscript, and scholarship boundaries

The package’s RNA and cross-species contracts distinguish executable readiness,
cohort readiness, descriptive analysis, and biological inference. Readable files
alone do not promote a cohort: current runtime/config/reference provenance,
manifest identity, successful quantification receipt, lock ownership, and
recovery state remain required. Positional cross-species alignment is not an
accepted fallback.

No manuscript claim, statistical result, figure, caption, citation ledger, or
rendered PDF was promoted by this review. The BeeWAS nested manuscript scaffold
must be reviewed and published independently, with its own current data,
provenance, citation verification, visual QA, and owner approval. A green
software suite does not establish statistical validity, accessibility conformance,
scientific completion, or publication readiness.

## Remaining gates and next actions

1. Recheck the hosted PR matrix after the latest parent push and investigate any
   newly failing lane from its fresh logs.
2. Run clean-wheel installation and Python 3.11/3.12 smoke tests in disposable
   environments; retain exact revisions and artifact hashes.
3. Reconcile and publish only intended BeeWAS changes through its nested PR;
   update the parent gitlink only after remote reachability is verified.
4. Reconcile the pre-existing Hymenoptera dirty snapshot on its own branch;
   do not stage or publish it from the parent.
5. Before any scientific completion statement, obtain a quiescent producer,
   current manifest, provenance-qualified finalized matrices, and a reproducible
   evidence bundle. Until then, biological inference remains withheld.
