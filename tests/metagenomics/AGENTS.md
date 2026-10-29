# AGENTS.md — `MetaInformAnt/tests/metagenomics/`

tests for the metagenomics module.
Files (verified 2026-10-28): __init__.py, test_metagenomics_amplicon.py, test_metagenomics_annotation_measures.py, test_metagenomics_bin_refinement.py, test_metagenomics_comparative.py, test_metagenomics_comparative_undefined.py, test_metagenomics_differential_abundance_mc.py, test_metagenomics_diversity.py, test_metagenomics_functional.py, test_metagenomics_input_checks.py, test_metagenomics_ordination.py, test_metagenomics_pathways_differential.py, test_metagenomics_scaffold_orientation.py, test_metagenomics_shotgun.py, test_metagenomics_student_t.py, test_metagenomics_visualization.py.

## Coverage by file
- `test_metagenomics_amplicon.py`: OTU clustering (the exact k-mer pre-filter, `META_OTU_THRESHOLD` resolution), chimera detection (`min_divergence`), ASV denoising (Poisson tail, Phred quality scores, the SciPy-free path), error-rate estimation, paired-read merging, naive Bayes taxonomy (abstention, seeded bootstrap) and taxonomy tree counts.
- `test_metagenomics_comparative.py`: differential abundance (direction, NumPy inputs), CLR transform, indicator species, effect sizes, biomarker discovery (`cv_folds` capped by the smallest class) and the SciPy-free Kruskal-Wallis and chi-square tail.
- `test_metagenomics_student_t.py`: the exact Student's t tail (regularized incomplete beta) against closed forms (df 1 and 2) and SciPy, and the SciPy-free Welch test (run in a child interpreter in which SciPy cannot be imported).
- `test_metagenomics_comparative_undefined.py`: undefined tests report NaN rather than p = 1.0 (Welch, Kruskal-Wallis), Benjamini-Hochberg over the tests that ran, the ANCOM-like W statistic (no p-value), `cv_accuracy` of `biomarker_discovery` being `None` with a reason, and the `aldex2_like` Monte Carlo without NumPy.
- `test_metagenomics_differential_abundance_mc.py`: the seeded Dirichlet Monte Carlo of `aldex2_like`.
- `test_metagenomics_diversity.py`: alpha and beta diversity (Fisher's alpha, integer-count validation, NumPy inputs), rarefaction, PERMANOVA, ordination and the seeded NMDS random starts.
- `test_metagenomics_functional.py`: ORF prediction (translation tables, each ORF reported once), annotation (no fabricated E-values or family labels) and pathway completeness by steps.
- `test_metagenomics_annotation_measures.py`: `profile_coverage` and `consensus_identity` of the best hit, real hit bounds, and that nothing to align is never a hit.
- `test_metagenomics_bin_refinement.py`: `refine_bins` applies the contamination criterion to kept bins and split sub-bins (GC content and coverage of sub-bins), and `bin_contigs` input checks.
- `test_metagenomics_ordination.py`: PCoA with negative eigenvalues (exact eigensolver, honest `variance_explained`), Kruskal's NMDS (stress, planted-configuration recovery after Procrustes alignment, stress against scikit-learn, seeding), `plot_ordination` labels and distance checks, and PERMANOVA on perfectly separated groups.
- `test_metagenomics_scaffold_orientation.py`: scaffold order, orientation and gaps recovered from simulated FR read pairs of a genome cut into shuffled, flipped contigs, and the reported limitation for contigs that are their own reverse complement.
- `test_metagenomics_input_checks.py`: Krona chart normalisation, strand-independent k-mer profiling, hash-seed-independent read merging, error-model input checks, the pathway database check and placeholder-free estimators.
- `test_metagenomics_pathways_differential.py`: `find_differential_pathways`.
- `test_metagenomics_shotgun.py`: deterministic assembly, scaffolding on both strands, TNF binning, marker-gene bin quality, high-quality bin counts, refinement with markers, and k-mer profiling when taxa share no ancestor.
- `test_metagenomics_visualization.py`: every plot function, and `plot_ordination` method selection (PCoA or NMDS).

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/metagenomics` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.
