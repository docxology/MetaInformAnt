# METAINFORMANT 0.4.0 Migration Guide

Version 0.4.0 is a deliberate pre-1.0 breaking release. Update callers before
using the package in a production or scientific workflow.

## Removed import paths

Use the canonical modules below. The former compatibility modules are removed:

| Removed | Replacement |
| --- | --- |
| `metainformant.core.config` | `metainformant.core.utils.config` |
| `metainformant.core.paths` | `metainformant.core.io.paths` |
| `metainformant.dna.sequences` | `metainformant.dna.sequence.core` |

## Removed command wrappers

Use `bash scripts/package/test.sh` and its documented `--mode` options. The
deprecated `scripts/package/run_tests.sh` and
`scripts/package/uv_test_optimized.sh` wrappers are removed.

## RNA and external data behavior

NCBI-facing operations require either an explicit contact email or
`NCBI_EMAIL`. Anonymous requests require an explicit `allow_anonymous=True`
opt-in. Anonymous mode emits a warning, does not invent or persist an email,
and records the contact mode in provenance.

RNA subprocesses use campaign-local `NCBI_SETTINGS`, `TMPDIR`, `VDB_CONFIG`,
and cache paths. The package no longer modifies user-global SRA Toolkit
configuration or scans home-directory/tmp SRA locations during cleanup.

## Cross-species analysis

Cross-species comparisons require matching labels or an explicit validated
alignment map. Positional fallback is removed. Comparisons with incomplete
metadata remain descriptive and cannot be promoted to biological inference.
