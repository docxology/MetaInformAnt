# dna tests

pytest suite for the `dna` domain of METAINFORMANT. Tests import from `src/metainformant/dna` and follow the repo's real-implementation policy.

## Files

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

Subdirectories: `data`.

Run from the repo root: `pytest tests/dna/ -v`. Tests follow the real-implementation policy (no mocks).
