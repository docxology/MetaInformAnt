# Phenotype Module

The `phenotype` module provides tools for phenotypic trait analysis, curation, and integration with genotypic data.

## Overview

This module handles morphological and behavioral phenotype data, including automated curation from web sources like AntWiki.

## Submodules

### AntWiki Integration (`antwiki.py`)
Automated phenotype data collection from AntWiki database.

**Key Features:**
- Species trait extraction and parsing
- Morphological measurement standardization
- Behavioral trait classification
- Data validation and quality control

**Usage:**
```python
from metainformant.phenotype import antwiki

# Retrieve phenotype data
species_traits = antwiki.get_species_traits("Camponotus pennsylvanicus")
morphological_data = antwiki.extract_morphological_measurements(species_traits)
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import population
from metainformant.phenotype import antwiki

# Genotype-phenotype association analysis
genetic_markers = population.get_polymorphic_sites(sequences)
phenotypic_traits = antwiki.get_phenotypic_variation(species)
associations = analyze_genotype_phenotype_associations(genetic_markers, phenotypic_traits)
```

### With Ontology Module
```python
from metainformant.phenotype import antwiki
from metainformant.ontology import go

# Functional annotation of phenotypic traits
traits = antwiki.get_behavioral_traits(species)
functional_annotation = go.analyze_trait_ontology(traits)
```

## Data Sources

- AntWiki database for ant species phenotypes
- Morphological measurement databases
- Behavioral observation datasets
- Quantitative trait databases

## Performance Features

- Efficient web scraping with rate limiting
- Caching of retrieved phenotype data
- Batch processing for multiple species
- Incremental data updates

## Testing

Comprehensive tests cover:
- Web scraping functionality
- Data parsing accuracy
- Trait standardization
- Integration with other modules

## Dependencies

- BeautifulSoup for HTML parsing
- Requests for HTTP operations
- Optional: species-specific databases

This module provides essential tools for phenotype-genotype association studies and morphological analysis.
