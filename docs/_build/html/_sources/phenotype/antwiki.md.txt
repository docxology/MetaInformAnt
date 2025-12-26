### Phenotype: AntWiki

#### Loading AntWiki JSON Data

Function: `load_antwiki_json`

```python
from pathlib import Path
from metainformant.phenotype import antwiki

entries = antwiki.load_antwiki_json(Path("tests/data/phenotype/antwiki_dataset_sorted_final_01.json"))
```

#### Web Scraping AntWiki

Comprehensive web scraping functionality for extracting species data from AntWiki.

**Basic Usage:**
```python
from metainformant.phenotype.scraper import AntWikiScraper, load_scraper_config

# Load configuration
config = load_scraper_config()

# Initialize scraper
scraper = AntWikiScraper(config)

# Scrape single species
data = scraper.scrape_species_page("Camponotus_pennsylvanicus")
print(f"Species: {data['species']}")
print(f"Traits: {data['traits']}")
print(f"Measurements: {data['measurements']}")
```

**Scraping All Species:**
```python
from pathlib import Path

# Scrape all species with progress tracking
stats = scraper.scrape_all_species(
    output_dir=Path("output/phenotype/antwiki/"),
    limit=10,  # Optional: limit for testing
    resume=True  # Resume from checkpoint
)

print(f"Completed: {stats['completed']}, Failed: {stats['failed']}")
```

**Command-line Interface:**
```bash
# Scrape single species
python3 scripts/phenotype/scrape_antwiki.py --species Camponotus_pennsylvanicus

# Scrape all species
python3 scripts/phenotype/scrape_antwiki.py --all --limit 10

# Custom delay and output directory
python3 scripts/phenotype/scrape_antwiki.py --all --delay 3.0 --output output/phenotype/antwiki/
```

**Extracted Data:**
- Measurements: Morphological measurements from tables and infoboxes
- Traits: Behavioral and ecological traits from multiple sections
- Taxonomy: Taxonomic hierarchy information
- Distribution: Geographic distribution data
- Description: Full species description from main content
- Images: Image URLs from the page

**Configuration:**
Configuration is loaded from `config/phenotype/antwiki_scraper.yaml`:
- Rate limiting delay (default: 2 seconds, configurable via `PHEN_SCRAPE_DELAY`)
- Robots.txt checking (enabled by default)
- Request timeout and retry settings
- Output directory defaults

**Rate Limiting:**
The scraper respects rate limiting with configurable delays between requests. Default is 2 seconds, which can be overridden via configuration or command-line arguments.
