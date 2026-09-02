# Agent Directives: src/metainformant/eqtl

Module `eqtl` of the METAINFORMANT package. Extracted from
`scripts/eqtl/` (thin-orchestrator migration, 2026-09-01).

## Role
Reusable eQTL input construction, transcriptome SNP-calling wrappers
(HISAT2/samtools/bcftools), bcftools stats parsing, and pipeline
orchestration. Scripts are thin; business logic lives here.

## Contents
- `synthetic.py` - synthetic data generators + real quantification loader
- `variant_calling.py` - external-tool wrappers (idempotent, resume-safe)
- `variant_stats.py` - bcftools stats parsing and summaries
- `pipeline.py` - parameter resolution + end-to-end run

## Rules
- Scripts under `scripts/eqtl/` must only orchestrate; add methods here.
- External tools are subprocess-wrapped; check availability with
  `check_tools()` before runs; never mock tools in tests.
- Keep real-expression loading separate from synthetic generation so
  provenance stays explicit (real vs synthetic genotypes).
- Repo-wide policy: see the repository-root `AGENTS.md`.
