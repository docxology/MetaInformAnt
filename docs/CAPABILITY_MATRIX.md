# MetaInformAnt capability matrix

This is the source-linked capability contract for the 0.4.0 release. A
function being importable does not make it a completed scientific workflow.
The matrix describes implementation readiness; biological inference requires a
separate current manifest, provenance-qualified inputs, statistical evidence,
and an evidence bundle.

| Surface | Source of truth | Capability | Release boundary |
| --- | --- | --- | --- |
| Core utilities, I/O, logging, configuration | `src/metainformant/core/` | Stable | Public package utilities and canonical imports only |
| DNA sequence and genome helpers | `src/metainformant/dna/` | Stable | NCBI access requires explicit contact or anonymous opt-in |
| RNA acquisition and Amalgkit orchestration | `src/metainformant/rna/engine/`, `src/metainformant/rna/amalgkit/` | Experimental / operational | Requires pinned Amalgkit capability, campaign-local SRA environment, locks, receipts, and resumable state |
| RNA quantification provenance and readiness | `src/metainformant/rna/engine/provenance.py`, `scripts/rna/check_pipeline_status.py` | Stable contract | Readable files without current provenance are not promoted |
| RNA cross-species comparison | `src/metainformant/rna/analysis/cross_species.py` | Experimental | Explicit sample alignment, orthology validation, and descriptive statistics are required |
| GWAS and BeeWAS client/reporting surfaces | `src/metainformant/gwas/`, `projects/apis_gwas/` | Experimental / nested release | Nested repository PR and provenance validators govern publication |
| General domain analysis modules | `src/metainformant/<module>/` | Stable or experimental per module README | API availability is not evidence of domain-complete analysis |
| Optional integrations and external services | Module-specific `README.md`, `SPEC.md`, and dependency declarations | Experimental | Missing optional dependencies must follow the documented contract |
| Scaffolds and placeholders | Source TODOs, explicit `NotImplementedError`, and module validation reports | Scaffold | Must not be presented as production or inference-ready |
| Removed 0.4.0 surfaces | `docs/MIGRATION_0.4.md` | Removed | No runtime compatibility shims; migration documentation is the sole reference |
| Descriptive matrices and permutation scores | RNA analysis outputs and their provenance | Descriptive only | Must not be labeled conserved/divergent biological inference without the required evidence |
| Biological inference and release claims | Current manifests, receipts, finalized matrices, and evidence bundle | Gated | Withheld while the producer is active or the cohort is partial |

## Promotion rules

1. Stable means the source, tests, documentation, and import contract agree;
   it does not guarantee external service availability.
2. Experimental means the implementation is usable under the documented
   contract but remains subject to fixture, external-tool, or scientific
   validation limits.
3. Scaffold means the surface is intentionally incomplete and cannot be used
   as evidence of a completed analysis.
4. A release claim may advance from executable readiness to cohort readiness,
   descriptive analysis, and biological inference only when each preceding
   state has its own current receipt and evidence.
