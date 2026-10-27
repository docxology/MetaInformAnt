# AGENTS.md — `MetaInformAnt/tests/phenotype/`

tests for the phenotype module.
Files (verified 2026-10-26): __init__.py, test_antwiki_records.py, test_antwiki_scraper_robots.py, test_mappings.py, test_multivariate.py, test_phenotype_basic.py, test_phenotype_behavior.py, test_phenotype_comprehensive.py, test_phenotype_integration.py, test_phenotype_life_course.py, test_phenotype_modules.py, test_phenotype_morphological.py, test_phenotype_scraper.py, test_phenotype_workflow.py, test_phewas_regression.py, test_phewas_small_sample_pvalues.py, test_plots.py, test_plots_subfigure_axes.py, test_scrape_antwiki_script.py, test_statistical.py (20 files).

## Coverage notes
- `test_antwiki_scraper_robots.py`: robots.txt rules loaded through `AntWikiScraper.load_robots_txt` (refusal before any request, per-origin scoping, user-agent groups, opt-out, redirect targets, `Crawl-delay` and retry backoff), config aliases and JSON config files, explicit species lists with `failed_species`, and output file names derived from page titles. Nothing opens a socket.
- `test_scrape_antwiki_script.py`: runs `scripts/phenotype/scrape_antwiki.py` as a subprocess on the paths that stop before any request (`--help`, `--all`, `--resume`, `--delay`).
- `test_phewas_small_sample_pvalues.py`: `run_phewas` p-values equal `scipy.stats.linregress` p-values at 10, 20 and 50 samples.
- `test_plots_subfigure_axes.py`: the phenotype plot functions return the root `Figure` (savable as PNG) for axes created through `fig.subfigures(...)` and their own figure for ordinary axes; categorical proportions are coloured from `tab20`.

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/phenotype -p no:cacheprovider` (verified 2026-09-05).
Repo-wide policy: see the repository-root `AGENTS.md`.