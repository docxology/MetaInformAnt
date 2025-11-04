"""Phenotype curation and scraping utilities."""

from .antwiki import load_antwiki_json

try:
    from .life_course import (
        aggregate_temporal_phenotypes,
        extract_phenotypes_from_events,
        map_events_to_traits,
    )

    __all__ = [
        "load_antwiki_json",
        "extract_phenotypes_from_events",
        "aggregate_temporal_phenotypes",
        "map_events_to_traits",
    ]
except ImportError:
    # life_course depends on life_events which may not be available
    __all__ = [
        "load_antwiki_json",
    ]
