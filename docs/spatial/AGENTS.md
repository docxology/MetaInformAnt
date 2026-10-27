# Agent Directives: docs/spatial

## Role

Maintain spatial transcriptomics documentation for tissue-coordinate expression
analysis.

## Scope

- Document current interfaces under `src/metainformant/spatial/`.
- Link to local spatial topic files such as [index.md](index.md), [io.md](io.md),
  [clustering.md](clustering.md), and [visualization.md](visualization.md).
- Keep examples runnable or clearly marked as pseudocode.

## Boundaries

- Use `metainformant.core.io` for domain data file I/O in examples.
- Direct stdlib parsing is acceptable only for protocol glue, CLI helpers, core
  utilities, or narrow parser internals covered by tests.
- Do not claim support for external spatial platforms unless source code or
  tests in this checkout exercise that path.
- Place generated figures and analysis output under `output/`.

## Documentation Rules

- Every call quoted in a `python` or `python-snippet` block must bind to the real signature; run the documentation gate (`scripts/verify_documentation_code.py`) after editing.
- Say which coordinate order a snippet assumes: Visium `[row, col]` (`"row_col"`), Xenium and MERFISH `[x, y]` (`"xy"`); plots of bare arrays need `coordinate_order="xy"` for the latter two.
- Keep `CAPABILITIES.md` entries in step with the docstrings and signatures of the functions they mirror.
- Name the seed behaviour of any function that defaults to `seed=None` (`compute_ligand_receptor_interactions`, `niche_identification`).
- Show the gene names whenever a snippet scores ligand-receptor pairs or checks deconvolution markers (`gene_names=...`, or a dataset or DataFrame that carries them): genes are matched by exact, case-sensitive name, and the built-in ligand-receptor database is written in human HGNC symbols.
