# AGENTS.md — `MetaInformAnt/tests/metagenomics/`

tests for the metagenomics module.
Files (verified 2026-10-27): __init__.py, test_metagenomics_amplicon.py, test_metagenomics_comparative.py, test_metagenomics_differential_abundance_mc.py, test_metagenomics_diversity.py, test_metagenomics_functional.py, test_metagenomics_pathways_differential.py, test_metagenomics_shotgun.py, test_metagenomics_visualization.py.

## Coverage by file
- `test_metagenomics_amplicon.py`: OTU clustering (the exact k-mer pre-filter, `META_OTU_THRESHOLD` resolution), chimera detection (`min_divergence`), ASV denoising (Poisson tail, Phred quality scores, the SciPy-free path), error-rate estimation, paired-read merging, naive Bayes taxonomy (abstention, seeded bootstrap) and taxonomy tree counts.
- `test_metagenomics_comparative.py`: differential abundance (direction, NumPy inputs), CLR transform, indicator species, effect sizes, biomarker discovery (`cv_folds` capped by the smallest class) and the SciPy-free Kruskal-Wallis and chi-square tail.
- `test_metagenomics_differential_abundance_mc.py`: the seeded Dirichlet Monte Carlo of `aldex2_like`.
- `test_metagenomics_diversity.py`: alpha and beta diversity (Fisher's alpha, integer-count validation, NumPy inputs), rarefaction, PERMANOVA, ordination and seeded NMDS.
- `test_metagenomics_functional.py`: ORF prediction (translation tables, each ORF reported once), annotation (no fabricated E-values or family labels) and pathway completeness by steps.
- `test_metagenomics_pathways_differential.py`: `find_differential_pathways`.
- `test_metagenomics_shotgun.py`: deterministic assembly, scaffolding on both strands, TNF binning, marker-gene bin quality, high-quality bin counts, refinement with markers, and k-mer profiling when taxa share no ancestor.
- `test_metagenomics_visualization.py`: every plot function, and `plot_ordination` method selection (PCoA or NMDS).

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/metagenomics` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.
