# Manuscript status — MetaInformAnt

## Repo type

Tool/library repository: METAINFORMANT is a comprehensive bioinformatics
toolkit for multi-omic analysis (per the root `README.md`: 28 analysis modules
under `src/metainformant/`, CLI- and API-facing, Apache-2.0 licensed).

## Why no publication-target manuscript applies today

Files checked in making this determination: root `README.md` (toolkit
positioning), `pyproject.toml` (package metadata), `docs/` (module
documentation, capability matrix, tutorials — a software documentation tree,
not a research write-up), and `SPEC.md`/`TODO.md` (engineering scope). There
are no research analysis outputs, figure registries, or claim ledgers in the
repo that a publication-track manuscript could describe. Per the
docxology/template standard, tool repositories get a status marker, not
fabricated section stubs.

Note: a research manuscript for the companion RNA-seq workflow lives in the
`hymenoptera_amalgkit` submodule
(`projects/hymenoptera_amalgkit/doc/manuscript/`, including
`reproducibility_checklist.md`); it is owned by that submodule, not this repo.
`docs/ORCHESTRATION.md` links into it and resolves only when the submodule is
initialized.

## What would trigger creating a manuscript

- A methods paper describing the toolkit's benchmarked analysis modules, with
  real benchmark outputs and figures committed under `output/`.
- An owner decision to write up a specific applied analysis performed with
  MetaInformAnt.

Then a full `manuscript/` (numbered sections, `config.yaml`, `preamble.md`,
`references.bib`) should be created at the repo top level.
