# SPEC: Phenotype Scripts

Data acquisition and analysis scripts for phenotypic traits and life course events.

## Pipelines

### 1. AntWiki Scraping
Automated retrieval and parsing of taxonomic and phenotypic data from AntWiki.
- `scrape_antwiki.py`: Orchestrates the crawl, parse, and JSON-save cycle.

### 2. Phenotype Extraction
Deriving phenotypic traits from life course event sequences.
- `extract_phenotypes.py`: Maps events to specific trait categories.

## Standards

- **Validation**: Scraped data must pass the schema validation defined in `core.validation`.
- **Anonymization**: All individual-level data must be anonymized or aggregated before export.
