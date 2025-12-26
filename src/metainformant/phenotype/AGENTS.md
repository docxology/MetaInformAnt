# AI Agents in Phenotype Development

This document outlines AI assistance in developing METAINFORMANT's phenotypic trait analysis capabilities.

## AI Contributions

### Phenotype Architecture
**Code Assistant Agent** designed:
- AntWiki JSON data loading with validation
- Life course phenotype extraction from event sequences
- Temporal phenotype aggregation utilities
- Integration with life_events module

### Core Components
**Code Assistant Agent** implemented:

#### AntWiki Integration (`antwiki.py`)
- `load_antwiki_json()`: Loads and validates AntWiki JSON format data
- Uses core.io utilities for robust file handling (gzip support, error handling)
- Data validation for expected structure (species/taxon, measurements, traits)
- Proper error handling with core.errors (IOError, ValidationError)
- Logging integration with core.logging

#### AntWiki Web Scraping (`scraper.py`)
- `AntWikiScraper`: Comprehensive web scraper for AntWiki species pages
- `load_scraper_config()`: Configuration loading with environment variable overrides
- `scrape_species_page()`: Scrape single species page with all sections
- `scrape_all_species()`: Batch scraping with progress tracking and checkpoint/resume
- `get_species_list()`: Discover all species pages from AntWiki
- Data extraction methods:
  - `extract_measurements()`: Parse morphological measurements from tables/infoboxes
  - `extract_traits()`: Extract behavioral/ecological traits from multiple sections
  - `extract_taxonomy()`: Extract taxonomic information
  - `extract_distribution()`: Extract geographic distribution data
  - `extract_description()`: Extract species description
  - `extract_images()`: Extract image URLs
- Rate limiting with configurable delay (default: 2 seconds)
- Robots.txt compliance checking
- Retry logic with exponential backoff using core.errors.retry_with_backoff
- Progress tracking and checkpoint/resume support
- Organized output structure (individual species files + combined dataset)
- Uses BeautifulSoup4 for HTML parsing
- Real HTTP requests with requests library (no mocks)

#### Life Course Analysis (`life_course.py`)
- `extract_phenotypes_from_events()`: Extracts phenotypes from single EventSequence
- `aggregate_temporal_phenotypes()`: Aggregates phenotypes across multiple sequences in time windows
- `map_events_to_traits()`: Maps events to predefined trait categories
- Integration with life_events module for temporal event analysis
- Input validation and error handling
- Logging for operations and debugging

### Quality Assurance
**Code Assistant Agent** contributed to:
- Comprehensive error handling with core utilities
- Input validation for all functions
- Logging for debugging and monitoring
- Type hints and documentation
- Integration with core.io, core.errors, core.logging

## Development Approach

- **Core Integration**: All functions use core utilities (io, logging, errors, paths)
- **Error Handling**: Proper exception hierarchy with ValidationError and IOError
- **Data Validation**: Validates data structure and input parameters
- **Optional Dependencies**: Graceful handling of optional life_events module
- **Documentation**: Comprehensive docstrings with examples

## Quality Assurance

- Human oversight ensures phenotype accuracy and biological relevance
- AI assistance accelerates development while maintaining standards
- Comprehensive testing validates phenotype functionality
- Real implementations without mocking (following project standards)

## Current Capabilities

### AntWiki Data Loading
- Load JSON files (supports both list and single dict formats)
- Validate data structure (species/taxon, measurements, traits)
- Handle errors gracefully (FileNotFoundError, IOError, ValidationError)
- Support gzip-compressed files via core.io

### AntWiki Web Scraping
- Comprehensive scraping of all species pages and all sections
- Extract measurements, traits, taxonomy, distribution, descriptions, and images
- Rate limiting with configurable delay (respects robots.txt)
- Retry logic with exponential backoff for network failures
- Progress tracking with checkpoint/resume support
- Organized output structure (individual species files + combined dataset)
- Command-line interface for easy usage
- Configuration via YAML with environment variable overrides

### Life Course Phenotype Extraction
- Extract phenotypes from event sequences (requires life_events module)
- Aggregate temporal patterns across multiple sequences
- Map events to trait categories
- Handle empty sequences and invalid inputs gracefully

This phenotype infrastructure provides a solid foundation for phenotypic trait analysis with proper error handling, validation, and integration with core utilities.

## Complete Function Signatures

### AntWiki Data Loading (`antwiki.py`)
- `load_antwiki_json(path: Path, validate: bool = True) -> list[dict[str, Any]]`

### AntWiki Web Scraping (`scraper.py`)
- `load_scraper_config(config_path: Path | None = None) -> AntWikiScraperConfig`

### Life Course Analysis (`life_course.py`)
- `extract_phenotypes_from_events(event_sequence: Any, trait_mapping: dict[str, list[str]] | None = None) -> dict[str, Any]`
- `aggregate_temporal_phenotypes(sequences: list[Any], time_windows: list[tuple[float, float]], trait_categories: list[str]) -> dict[str, Any]`
- `map_events_to_traits(event_sequence: Any, trait_definitions: dict[str, list[str]]) -> dict[str, list[Any]]`
