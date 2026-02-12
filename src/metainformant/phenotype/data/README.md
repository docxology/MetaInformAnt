# Phenotype Data

AntWiki data loading, validation, scraping, and phenotype extraction for myrmecological trait analysis.

## Contents

| File | Purpose |
|------|---------|
| `antwiki.py` | AntWiki record loading, filtering, phenotype matrices, and reporting |
| `scraper.py` | Web scraper for AntWiki species pages with configurable extraction |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `AntWikiRecord` | Dataclass for a species record with taxonomy and phenotype fields |
| `load_antwiki_json()` | Load and validate AntWiki JSON data files |
| `save_antwiki_json()` | Serialize AntWikiRecord list to JSON |
| `filter_antwiki_records()` | Filter records by taxonomy, region, or trait criteria |
| `get_phenotype_distribution()` | Distribution statistics for a named phenotype |
| `find_similar_species()` | Find species with similar phenotype profiles |
| `create_phenotype_matrix()` | Species-by-trait matrix for multivariate analysis |
| `AntWikiScraper` | Configurable scraper for AntWiki species pages |
| `AntWikiScraperConfig` | Configuration dataclass for scraping parameters |

## Usage

```python
from metainformant.phenotype.data.antwiki import load_antwiki_json, create_phenotype_matrix
from metainformant.phenotype.data.scraper import AntWikiScraper

records = load_antwiki_json("data/antwiki_records.json")
matrix = create_phenotype_matrix(records)
scraper = AntWikiScraper(output_dir="output/phenotype")
```
