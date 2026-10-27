# AGENTS.md — `MetaInformAnt/tests/phenotype/`

tests for the phenotype module.
Files (verified 2026-10-27): __init__.py, test_antwiki_records.py, test_antwiki_scraper_robots.py, test_cross_omic_pvalues.py, test_cross_omic_tied_ranks.py, test_life_course_no_events.py, test_life_course_outcome_analysis.py, test_mappings.py, test_morphological_unit_conversion.py, test_multivariate.py, test_phenotype_basic.py, test_phenotype_behavior.py, test_phenotype_comprehensive.py, test_phenotype_integration.py, test_phenotype_life_course.py, test_phenotype_modules.py, test_phenotype_morphological.py, test_phenotype_scraper.py, test_phenotype_workflow.py, test_phewas_regression.py, test_phewas_small_sample_pvalues.py, test_plots.py, test_plots_subfigure_axes.py, test_scrape_antwiki_script.py, test_statistical.py (24 test files).

## Coverage notes
- `test_antwiki_scraper_robots.py`: robots.txt rules loaded through `AntWikiScraper.load_robots_txt` (refusal before any request, per-origin scoping, user-agent groups, opt-out, redirect targets, `Crawl-delay` and retry backoff), config aliases and JSON config files, explicit species lists with `failed_species`, and output file names derived from page titles. Nothing opens a socket.
- `test_scrape_antwiki_script.py`: runs `scripts/phenotype/scrape_antwiki.py` as a subprocess on the paths that stop before any request (`--help`, `--all`, `--resume`, `--delay`).
- `test_phewas_small_sample_pvalues.py`: `run_phewas` p-values equal `scipy.stats.linregress` p-values at 10, 20 and 50 samples.
- `test_plots_subfigure_axes.py`: the phenotype plot functions return the root `Figure` (savable as PNG) for axes created through `fig.subfigures(...)` and their own figure for ordinary axes; categorical proportions are coloured from `tab20`.

## Contract tests

- `test_cross_omic_pvalues.py`: the p-values of `cross_omic` come from the exact Student's t and F distributions (`scipy.stats.t.sf` and `scipy.stats.f.sf`), at any sample size.
- `test_cross_omic_tied_ranks.py`: Spearman correlation averages tied ranks, equals `scipy.stats.spearmanr` on data with ties and does not depend on the order of the samples.
- `test_life_course_outcome_analysis.py`: `analyze_life_course` treats a sequence without `metadata["outcomes"]` as unknown, reports counts and prevalence among sequences with data, compares event count, duration and complexity with Mann-Whitney, and tests each event type with Fisher's exact test.
- `test_life_course_no_events.py`: `n_sequences_with_events` and `mean_duration` replace `analysis_timestamp`; with no events `mean_duration`, `duration_stats` and `complexity_stats` are `None` and `event_frequency` is `{}` instead of raising.
- `test_morphological_unit_conversion.py`: `Measurement.convert` is exact for every pair of `um` (also `µm` and `μm`), `mm`, `cm` and `m`, converts the uncertainty too, and raises `ValueError` for any other unit.

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/phenotype -p no:cacheprovider`
Repo-wide policy: see the repository-root `AGENTS.md`.
