# Phenotype Data

AntWiki data loading, validation, scraping, and phenotype extraction for myrmecological trait analysis.

## Contents

| File | Purpose |
|------|---------|
| `antwiki.py` | AntWiki record loading, filtering, phenotype matrices, and reporting |
| `scraper.py` | Web scraper for AntWiki species pages (`requests`-based, robots.txt-aware) with configurable extraction |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `AntWikiRecord` | Species record with taxonomy and phenotype fields, built from a dictionary that has `species_name` and `genus` |
| `load_antwiki_json()` | Load AntWiki JSON files as a list of dictionaries, validating their structure unless `validate=False` |
| `save_antwiki_json()` | Serialize AntWikiRecord list to JSON |
| `filter_antwiki_records()` | Filter records by genus, subfamily, minimum confidence, or required phenotypes |
| `get_phenotype_distribution()` | Distribution statistics for a named phenotype |
| `find_similar_species()` | Find species with similar phenotype profiles |
| `create_phenotype_matrix()` | Species-by-trait matrix for multivariate analysis |
| `AntWikiScraper` | Scraper for AntWiki species pages; consults robots.txt by default |
| `AntWikiScraperConfig` | Configuration class for scraping parameters (delay, timeout, retries, user agent, `check_robots`, `output_dir`) |
| `load_scraper_config()` | Read scraper settings from a JSON file (the defaults when no path is given) |
| `RobotsDisallowedError` | Raised when robots.txt forbids a URL or cannot be retrieved |

## Usage

```python
from metainformant.phenotype.data.antwiki import AntWikiRecord, create_phenotype_matrix, load_antwiki_json
from metainformant.phenotype.data.scraper import AntWikiScraper, AntWikiScraperConfig

entries = load_antwiki_json("data/antwiki_records.json")  # a list of dictionaries
records = [AntWikiRecord(entry) for entry in entries]  # each entry needs species_name and genus
species_names, phenotype_names, matrix = create_phenotype_matrix(records)
scraper = AntWikiScraper(AntWikiScraperConfig(output_dir="output/phenotype"))
```

`AntWikiScraper` takes an `AntWikiScraperConfig`; the output directory is a setting of that config.
