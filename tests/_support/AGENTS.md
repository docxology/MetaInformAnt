# tests/_support - Shared synthetic-data factories

## Purpose

Deterministic, typed builders for the synthetic data used across test files:
expression frames/matrices, ortholog tables, species trees, phenotype and
genotype maps, GWAS association rows, spatial cell layouts, and small real
files (metadata TSV, FASTA, config YAML). Import as:

    from tests._support.synth import make_expression_frame, ...

## Standards

- Factories are deterministic (seed-controlled) and return real
  numpy/pandas objects or write real files under caller-provided paths.
- Zero-mocks policy applies: factories never fake the unit under test.
- New shared builders belong here, not copy-pasted into tests/*/*.py.
- Refactors of existing tests should route their local ``_make_*``
  helpers through these factories (thin local wrappers are fine to
  preserve call-site signatures).

## Files

- ``synth.py`` - the factory library (single module, no subpackage sprawl)

## See also

- ``../AGENTS.md`` - test-suite standards
- ``scripts/test/run_module_tests.py`` - per-module test runner
