"""Phenotype curation and analysis utilities.

This module provides tools for phenotypic trait analysis, including:
- AntWiki JSON data loading with validation
- AntWiki web scraping for comprehensive species data
- Life course phenotype extraction from event sequences (requires life_events module)
- Temporal phenotype aggregation
- Integration with genotypic data for association studies

The life_course functions require the optional life_events module. If not available,
only the load_antwiki_json and scraping functions will be exported. The module handles this
gracefully with conditional imports.
"""

from .antwiki import load_antwiki_json
from .scraper import AntWikiScraper, AntWikiScraperConfig, load_scraper_config

try:
    from .life_course import (
        aggregate_temporal_phenotypes,
        extract_phenotypes_from_events,
        map_events_to_traits,
    )

    __all__ = [
        "load_antwiki_json",
        "AntWikiScraper",
        "AntWikiScraperConfig",
        "load_scraper_config",
        "extract_phenotypes_from_events",
        "aggregate_temporal_phenotypes",
        "map_events_to_traits",
    ]
except ImportError:
    # life_course depends on life_events which may not be available
    __all__ = [
        "load_antwiki_json",
        "AntWikiScraper",
        "AntWikiScraperConfig",
        "load_scraper_config",
    ]
