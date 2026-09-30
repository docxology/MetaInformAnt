# AGENTS.md — `MetaInformAnt/tests/dna`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `dna` domain of METAINFORMANT. Tests import from `src/metainformant/dna` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_dna_accession.py`
- `test_dna_alignment.py`
- `test_dna_codon_usage.py`
- `test_dna_compatibility_facades.py`
- `test_dna_comprehensive.py`
- `test_dna_consensus.py`
- `test_dna_distances.py`
- `test_dna_entrez_integration.py`
- `test_dna_fastq.py`
- `test_dna_gc_skew_tm.py`
- `test_dna_genomes.py`
- `test_dna_kmer.py`
- `test_dna_kmer_distances.py`
- … (+24 more test modules)
- Subdirectories: `data`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/dna/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
