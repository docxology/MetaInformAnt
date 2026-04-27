# Personal AI Infrastructure (PAI) - dna

## Context & Intent
- **Path**: `docs/dna`
- **Purpose**: Comprehensive DNA sequence analysis documentation covering sequence operations, alignment algorithms, phylogenetics, population genetics, variant calling, and genome retrieval.
- **Domain**: docs

## Virtual Hierarchy
- **Type**: Documentation
- **Parent**: `docs`

## Maintenance Notes
- **System**: Part of the METAINFORMANT Domain layer (DNA analysis).
- **Style**: Strict type hinting, no mocks in tests, reproducible examples.
- **Stability**: Core API stable; new algorithms added through submodules.

## AI Workflows
- **Modification**: Run functional tests in `tests/test_dna_*.py` before committing.
- **Documentation**: Update `SPEC.md` if architectural patterns change; add examples to `README.md`.

## Cross-References
- **Main docs**: [index.md](../dna/index.md) — Module overview and quick start
- **Sequence ops**: [sequence.md](sequences.md) — FASTA parser, reverse complement, translation
- **Alignment**: [alignment.md](../dna/alignment.md) — Needleman-Wunsch, Smith-Waterman, MUSCLE
- **Population genetics**: [population.md](../dna/population.md) — Tajima's D, Fst, π, LD
- **Applications**: [variant_calling.md](../dna/variants.md) — GATK, FreeBayes, VEP annotation
- **Related modules**: [RNA](../rna/index.md) for expression; [GWAS](../gwas/index.md) for association; [Quality](../quality/index.md) for FASTQ/BAM QC
