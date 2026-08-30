# AGENTS.md — `MetaInformAnt/tests/dna/`

tests for DNA/sequence/alignment/phylogeny modules.
Files (verified 2026-08-29): __init__.py, test_dna_accession.py, test_dna_alignment.py, test_dna_codon_usage.py, test_dna_compatibility_facades.py, test_dna_comprehensive.py, test_dna_consensus.py, test_dna_distances.py … (38 files).


## Children (documented on disk)
- data/

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/dna` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.