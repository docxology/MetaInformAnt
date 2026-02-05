"""Data loading and scraping utilities for phenotype data."""

from .antwiki import (
    AntWikiRecord,
    load_antwiki_json,
    save_antwiki_json,
    filter_antwiki_records,
    get_phenotype_distribution,
    find_similar_species,
    create_phenotype_matrix,
    generate_antwiki_report,
)
from .scraper import AntWikiScraper, AntWikiScraperConfig

__all__ = [
    "AntWikiRecord",
    "load_antwiki_json",
    "save_antwiki_json",
    "filter_antwiki_records",
    "get_phenotype_distribution",
    "find_similar_species",
    "create_phenotype_matrix",
    "generate_antwiki_report",
    "AntWikiScraper",
    "AntWikiScraperConfig",
]
