# Phenotype Module

The `phenotype` module provides tools for phenotypic trait analysis, curation, and integration with genotypic data.

## Overview

This module handles morphological and behavioral phenotype data, including automated curation from web sources like AntWiki. Enables phenotype-genotype association studies and morphological trait analysis.

## Key Components

### AntWiki Integration (`antwiki.py`)
Load phenotype data from AntWiki JSON files.

**Usage:**
```python
from metainformant.phenotype.antwiki import load_antwiki_json
from pathlib import Path

# Load AntWiki phenotype data
data = load_antwiki_json(Path("data/antwiki_species.json"))

# Each entry contains species phenotype information
for species_data in data:
    species_name = species_data.get("species", "unknown")
    measurements = species_data.get("measurements", {})
    traits = species_data.get("traits", [])
    
    print(f"{species_name}: {len(traits)} traits")
```

### Life Course Integration (`life_course.py`)
Extract and analyze temporal phenotypes from life event sequences.

**Usage:**
```python
from metainformant.phenotype import (
    extract_phenotypes_from_events,
    aggregate_temporal_phenotypes,
    map_events_to_traits
)

# Extract phenotypes from life event sequences
phenotypes = extract_phenotypes_from_events(event_sequences)

# Aggregate temporal phenotypes
aggregated = aggregate_temporal_phenotypes(phenotypes)

# Map events to trait categories
trait_mapping = map_events_to_traits(event_sequences)
```

**Integration with Life Events Module:**
```python
from metainformant.life_events import EventSequence, load_sequences_from_json
from metainformant.phenotype import extract_phenotypes_from_events

# Load life event sequences
sequences = load_sequences_from_json("life_events.json")

# Extract phenotypic traits from events
phenotypes = extract_phenotypes_from_events(sequences)
```

**Data Structure:**
AntWiki JSON files contain species entries with:
- `species`: Species name
- `measurements`: Morphological measurements (e.g., worker length, head width)
- `traits`: Behavioral and morphological trait classifications

**Example:**
```python
# Typical AntWiki JSON structure
[
    {
        "species": "Camponotus pennsylvanicus",
        "measurements": {
            "worker_length_mm": [6.0, 13.0],
            "head_width_mm": [1.8, 3.2]
        },
        "traits": ["arboreal", "carnivorous", "polygynous"]
    },
    ...
]
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import population
from metainformant.phenotype import load_antwiki_json

# Genotype-phenotype association analysis

# Analyze population genetics
diversity = population.nucleotide_diversity(sequences)

# Load phenotype data
phenotype_data = load_antwiki_json(Path("antwiki_species.json"))
# Extract traits for analysis
# See genotype-phenotype association analysis tools in other modules
```

### With Life Events Module
```python
from metainformant.life_events import load_sequences_from_json
from metainformant.phenotype import extract_phenotypes_from_events

# Extract phenotypes from temporal event sequences
sequences = load_sequences_from_json("life_events.json")
phenotypes = extract_phenotypes_from_events(sequences)
# Analyze temporal phenotype patterns
```

### With Ontology Module
```python
from metainformant.phenotype import load_antwiki_json
from metainformant.ontology import load_go_obo

# Functional annotation of phenotypic traits

# Load phenotype data
phenotype_data = load_antwiki_json(Path("antwiki_species.json"))

# Load GO for functional analysis
go_onto = load_go_obo("go-basic.obo")
# Use GO for trait functional annotation
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
