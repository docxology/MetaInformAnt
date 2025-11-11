# Phenotype Module

The `phenotype` module provides tools for phenotypic trait analysis, curation, and integration with genotypic data.

## Overview

This module handles morphological and behavioral phenotype data, including loading from structured JSON files (e.g., AntWiki format). Enables phenotype-genotype association studies and morphological trait analysis.

### Module Architecture

```mermaid
graph TB
    subgraph "Phenotype Module"
        AntWiki[antwiki<br/>AntWiki Integration]
        LifeCourse[life_course<br/>Life Course Analysis]
    end
    
    subgraph "Input Data"
        AntWikiJSON[AntWiki JSON]
        EventSeqs[Event Sequences]
    end
    
    subgraph "Other Modules"
        LifeEvents[life_events]
        GWAS_Mod[gwas]
        Networks_Mod[networks]
    end
    
    AntWikiJSON --> AntWiki
    EventSeqs --> LifeCourse
    AntWiki --> LifeCourse
    LifeCourse --> LifeEvents
    LifeCourse --> GWAS_Mod
    LifeCourse --> Networks_Mod
```

### Phenotype Analysis Workflow

```mermaid
flowchart LR
    Start[Input Data] --> Source{Data Source?}
    Source -->|AntWiki| LoadAntWiki[Load AntWiki JSON]
    Source -->|Events| LoadEvents[Load Event Sequences]
    LoadAntWiki --> Extract[Extract Phenotypes]
    LoadEvents --> Extract
    Extract --> Aggregate[Aggregate Traits]
    Aggregate --> Analyze[Analysis]
    Analyze --> Associate[GWAS Association]
    Associate --> Output[Results]
```

## Key Components

### AntWiki Integration (`antwiki.py`)
Load phenotype data from AntWiki JSON files.

**Usage:**
```python
from metainformant.phenotype.antwiki import load_antwiki_json
from metainformant.core.io import write_json
from metainformant.core.paths import expand_and_resolve
from pathlib import Path

# Load AntWiki phenotype data with error handling
try:
    data_path = expand_and_resolve("data/antwiki_species.json")
    data = load_antwiki_json(data_path)
    
    # Each entry contains species phenotype information
    for species_data in data:
        species_name = species_data.get("species") or species_data.get("taxon", "unknown")
        measurements = species_data.get("measurements", {})
        traits = species_data.get("traits", [])
        
        print(f"{species_name}: {len(traits)} traits, {len(measurements)} measurements")
        
except FileNotFoundError as e:
    print(f"File not found: {e}")
except Exception as e:
    print(f"Error loading data: {e}")
```

**Error Handling:**
The function raises `FileNotFoundError` for missing files, `IOError` from `core.errors` for file read issues, and `ValidationError` for invalid data structures.

### Life Course Integration (`life_course.py`)
Extract and analyze temporal phenotypes from life event sequences.

**Usage:**
```python
from metainformant.phenotype import (
    extract_phenotypes_from_events,
    aggregate_temporal_phenotypes,
    map_events_to_traits
)
from metainformant.life_events import EventSequence, Event, load_sequences_from_json
from datetime import datetime

# Create or load event sequences
# Note: extract_phenotypes_from_events takes a SINGLE EventSequence, not a list
sequence = EventSequence(
    person_id="person_001",
    events=[
        Event("diabetes", datetime(2020, 1, 1), "health"),
        Event("bachelors", datetime(2010, 6, 1), "education"),
    ]
)

# Extract phenotypes from a single event sequence
phenotypes = extract_phenotypes_from_events(sequence)
print(f"Total events: {phenotypes['total_events']}")
print(f"Domains: {phenotypes['domains']}")

# Aggregate temporal phenotypes from MULTIPLE sequences
sequences_list = [sequence, ...]  # List of EventSequence objects
aggregated = aggregate_temporal_phenotypes(sequences_list, time_window_years=5.0)
print(f"Total people: {aggregated['aggregates']['total_people']}")

# Map events to trait categories for a SINGLE sequence
trait_mapping = map_events_to_traits(sequence)
print(f"Health issues: {trait_mapping['health_issues']['count']}")
```

**Integration with Life Events Module:**
```python
from metainformant.life_events import EventSequence, load_sequences_from_json
from metainformant.phenotype import extract_phenotypes_from_events
from metainformant.core.paths import expand_and_resolve
from pathlib import Path

# Load life event sequences (returns list of EventSequence objects)
sequences = load_sequences_from_json(Path("data/life_events.json"))

# Extract phenotypic traits from each sequence
# Note: extract_phenotypes_from_events takes a SINGLE EventSequence
all_phenotypes = []
for sequence in sequences:
    try:
        phenotypes = extract_phenotypes_from_events(sequence)
        all_phenotypes.append(phenotypes)
    except Exception as e:
        print(f"Error processing sequence {sequence.person_id}: {e}")

print(f"Processed {len(all_phenotypes)} sequences")
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
from metainformant.phenotype.antwiki import load_antwiki_json
from metainformant.core.paths import expand_and_resolve
from pathlib import Path

# Genotype-phenotype association analysis

# Analyze population genetics
diversity = population.nucleotide_diversity(sequences)

# Load phenotype data with proper error handling
try:
    phenotype_path = expand_and_resolve("data/antwiki_species.json")
    phenotype_data = load_antwiki_json(phenotype_path)
    
    # Extract traits for analysis
    for entry in phenotype_data:
        species = entry.get("species") or entry.get("taxon")
        traits = entry.get("traits", [])
        # Use with genotype data for association analysis
        # See genotype-phenotype association analysis tools in other modules
except FileNotFoundError:
    print("Phenotype data file not found")
except Exception as e:
    print(f"Error loading phenotype data: {e}")
```

### With Life Events Module
```python
from metainformant.life_events import load_sequences_from_json, EventSequence
from metainformant.phenotype import extract_phenotypes_from_events, aggregate_temporal_phenotypes
from metainformant.core.paths import expand_and_resolve
from pathlib import Path

# Extract phenotypes from temporal event sequences
sequences = load_sequences_from_json(Path("data/life_events.json"))

# Extract phenotypes from each sequence (takes single EventSequence)
phenotypes_list = []
for sequence in sequences:
    phenotypes = extract_phenotypes_from_events(sequence)
    phenotypes_list.append(phenotypes)

# Aggregate temporal patterns across all sequences
aggregated = aggregate_temporal_phenotypes(sequences, time_window_years=5.0)
# Analyze temporal phenotype patterns
```

### With Ontology Module
```python
from metainformant.phenotype.antwiki import load_antwiki_json
from metainformant.ontology import load_go_obo
from metainformant.core.paths import expand_and_resolve
from pathlib import Path

# Functional annotation of phenotypic traits

# Load phenotype data
try:
    phenotype_path = expand_and_resolve("data/antwiki_species.json")
    phenotype_data = load_antwiki_json(phenotype_path)
    
    # Load GO for functional analysis
    go_path = expand_and_resolve("data/go-basic.obo")
    go_onto = load_go_obo(go_path)
    
    # Use GO for trait functional annotation
    for entry in phenotype_data:
        traits = entry.get("traits", [])
        # Map traits to GO terms for functional analysis
except Exception as e:
    print(f"Error in integration: {e}")
```

## Data Sources

- AntWiki JSON format for ant species phenotypes
- Morphological measurement data in structured JSON
- Behavioral observation datasets in JSON format
- Quantitative trait databases

## Error Handling

All functions use proper error handling with core utilities:

```python
from metainformant.phenotype.antwiki import load_antwiki_json
from metainformant.core.errors import IOError, ValidationError
from pathlib import Path

try:
    data = load_antwiki_json(Path("data.json"))
except FileNotFoundError:
    print("File not found")
except ValidationError as e:
    print(f"Invalid data format: {e}")
except IOError as e:
    print(f"I/O error: {e}")
```

## Data Validation

The `load_antwiki_json` function validates data structure by default:
- Ensures entries have `species` or `taxon` field
- Validates `measurements` is a dictionary if present
- Validates `traits` is a list if present

Disable validation by passing `validate=False`:

```python
data = load_antwiki_json(path, validate=False)
```

## Testing

Comprehensive tests cover:
- JSON loading and parsing accuracy
- Data structure validation
- Error handling for missing files and invalid data
- Integration with life_events module
- Integration with other modules

## Dependencies

- `metainformant.core` - Core utilities (io, logging, errors, paths)
- `metainformant.life_events` - Optional dependency for life course functions

This module provides essential tools for phenotype-genotype association studies and morphological analysis.
