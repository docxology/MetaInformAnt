# Phenotype Documentation

This directory contains comprehensive documentation for METAINFORMANT's phenotypic trait analysis and curation capabilities.

**⚠️ Package Management**: All setup and dependency installation uses `uv`. Use `uv venv`, `uv pip install`, `uv run`, `uv sync`, `uv add`, and `uv remove` - never use `pip` directly. See **[UV Setup Guide](../UV_SETUP.md)** for details.

## Overview

behavioral phenotype data analysis, including automated curation from web sources like AntWiki.

### Module Architecture

```mermaid
    subgraph "Phenotype Module"
        AntWikiantwikiAntwikiIntegration[antwiki_AntWiki Integration]
        LifeCourselifeCourseLifeCourseAnalysis[life_course_Life Course Analysis]
        BehaviorbehaviorBehavioralAnalysis[behavior_Behavioral Analysis]
        ChemicalchemicalChemotypes[chemical_Chemotypes]
        SonicsonicAcousticSignals[sonic_Acoustic Signals]
        MorphmorphologicalMorphometrics[morphological_Morphometrics]
        ElectronicelectronicTrackingData[electronic_Tracking Data]
    end
```
## Documentation Files

### Core Phenotype Tools
- **`index.md`**: Phenotype domain overview and module index
- **`antwiki.md`**: AntWiki database integration and phenotype curation

## Related Source Code

- See `src/metainformant/phenotype/` for implementation details
- See `tests/test_phenotype_*.py` for comprehensive test coverage
- See `src/metainformant/phenotype/README.md` for module-specific documentation

## Usage Examples

The phenotype domain supports trait analysis:

```python
from metainformant.phenotype import antwiki
from metainformant.phenotype.scraper import AntWikiScraper, load_scraper_config

# Load existing AntWiki JSON data
entries = antwiki.load_antwiki_json(Path("data/antwiki_species.json"))

# Or scrape AntWiki directly
config = load_scraper_config()
scraper = AntWikiScraper(config)
data = scraper.scrape_species_page("Camponotus_pennsylvanicus")
```

## Integration

Phenotype analysis integrates with:
- **DNA analysis** for genotype-phenotype associations
- **Ontology** for trait functional annotation
- **Statistical methods** for trait correlation analysis
- **Visualization** for trait distribution plotting

## Testing

Comprehensive tests ensure phenotype reliability:
- Web scraping functionality validation
- Data parsing and standardization
- Trait measurement accuracy
- Integration with other modules

## Contributing

When adding new phenotype functionality:
1. Update trait analysis documentation
2. Add comprehensive data validation tests
3. Ensure compatibility with phenotype databases
4. Update integration examples

This documentation provides complete coverage of METAINFORMANT's phenotype analysis capabilities.
