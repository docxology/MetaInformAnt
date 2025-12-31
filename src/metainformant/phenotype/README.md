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

### AntWiki Data Processing Pipeline

```mermaid
graph TD
    A[AntWiki Data Sources] --> B{Data Type}
    B -->|JSON Files| C[Load JSON Files]
    B -->|Web Scraping| D[Scrape AntWiki]

    C --> E[Parse JSON Structure]
    D --> F[HTML Parsing]

    E --> G[Extract Fields]
    F --> G

    G --> H[Validate Data]
    H --> I{Clean Data?}

    I -->|Yes| J[Data Cleaning]
    I -->|No| K[Raw Data]

    J --> L[Standardized Format]
    K --> L

    L --> M[Phenotype Categories]
    M --> N[Morphological Traits]
    M --> O[Behavioral Traits]
    M --> P[Ecological Traits]

    N --> Q[Measurements DB]
    O --> Q
    P --> Q

    Q --> R[Analysis Ready]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style L fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style R fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Data Fields"
        S[Species Info] -.-> G
        T[Measurements] -.-> G
        U[Traits] -.-> G
        V[Taxonomy] -.-> G
        W[Distribution] -.-> G
    end

    subgraph "Quality Control"
        X[Missing Values] -.-> H
        Y[Outliers] -.-> H
        Z[Data Types] -.-> H
        AA[Consistency] -.-> H
    end

    subgraph "Trait Types"
        BB[Morphometric] -.-> N
        CC[Coloration] -.-> N
        DD[Behavior] -.-> O
        EE[Ecology] -.-> P
    end
```

### Life Course Phenotype Analysis

```mermaid
graph TD
    A[Event Sequences] --> B[Event Parsing]
    B --> C[Extract Phenotypes]

    C --> D[Temporal Windows]
    D --> E[Window Aggregation]

    E --> F[Trait Categories]
    F --> G[Morphological Changes]
    F --> H[Behavioral Patterns]
    F --> I[Developmental Stages]

    G --> J[Growth Trajectories]
    H --> K[Activity Patterns]
    I --> L[Life Stage Transitions]

    J --> M[Statistical Analysis]
    K --> M
    L --> M

    M --> N[Phenotype Networks]
    N --> O[Association Mining]

    O --> P[Key Phenotypes]
    P --> Q[Biological Insights]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style E fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style Q fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Event Types"
        R[Developmental] -.-> B
        S[Morphological] -.-> B
        T[Behavioral] -.-> B
        U[Environmental] -.-> B
    end

    subgraph "Time Windows"
        V[Early Development] -.-> D
        W[Juvenile] -.-> D
        X[Adult] -.-> D
        Y[Senescence] -.-> D
    end

    subgraph "Analysis Methods"
        Z[Trajectory Modeling] -.-> M
        AA[Pattern Recognition] -.-> M
        BB[Network Analysis] -.-> N
        CC[Association Rules] -.-> O
    end
```

### Phenotype-Genotype Association Analysis

```mermaid
graph TD
    A[Phenotype Data] --> B[Genotype Data]
    B --> C[Sample Matching]

    A --> D[Phenotype Processing]
    D --> E[Normalize Traits]
    E --> F[Quality Control]

    C --> G[Association Testing]
    F --> G

    G --> H{Association Method}
    H -->|Linear Regression| I[Quantitative Traits]
    H -->|Logistic Regression| J[Binary Traits]
    H -->|ANOVA| K[Categorical Traits]
    H -->|Correlation| L[Continuous Traits]

    I --> M[Statistical Results]
    J --> M
    K --> M
    L --> M

    M --> N[Multiple Testing Correction]
    N --> O[Significant Associations]

    O --> P[Manhattan Plot]
    O --> Q[QQ Plot]
    O --> R[Regional Plot]

    P --> S[Genome-wide View]
    Q --> S
    R --> S

    S --> T[Candidate Genes]
    T --> U[Functional Validation]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style G fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style U fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Phenotype Types"
        V[Morphometric] -.-> D
        W[Behavioral] -.-> D
        X[Physiological] -.-> D
        Y[Developmental] -.-> D
    end

    subgraph "Genotype Data"
        Z[VCF Files] -.-> B
        AA[PLINK Format] -.-> B
        BB[Genotype Matrix] -.-> B
    end

    subgraph "Correction Methods"
        CC[Bonferroni] -.-> N
        DD[FDR] -.-> N
        EE[Permutation] -.-> N
    end
```

### Web Scraping and Data Collection

```mermaid
graph TD
    A[AntWiki Website] --> B[Discover Species]
    B --> C[Species List]

    C --> D{Scraping Strategy}
    D -->|Single Species| E[Scrape One Page]
    D -->|Batch Scraping| F[Scrape Multiple]

    E --> G[Parse HTML]
    F --> H[Queue Management]

    G --> I[Extract Content]
    H --> I

    I --> J[Data Sections]
    J --> K[Morphological Data]
    J --> L[Behavioral Data]
    J --> M[Ecological Data]
    J --> N[Taxonomic Data]

    K --> O[Structured Format]
    L --> O
    M --> O
    N --> O

    O --> P[Quality Validation]
    P --> Q{Save Results}
    Q -->|Individual Files| R[Per Species File]
    Q -->|Combined| S[Master Dataset]

    R --> T[Storage]
    S --> T

    T --> U[Analysis Pipeline]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style I fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style U fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Scraping Controls"
        V[Rate Limiting] -.-> D
        W[Error Handling] -.-> D
        X[Progress Tracking] -.-> F
        Y[Checkpoint/Resume] -.-> F
    end

    subgraph "Content Sections"
        Z[Measurements Table] -.-> K
        AA[Trait Descriptions] -.-> L
        BB[Distribution Maps] -.-> M
        CC[Taxonomy Info] -.-> N
    end

    subgraph "Output Formats"
        DD[JSON] -.-> R
        EE[CSV] -.-> S
        FF[Database] -.-> T
    end
```

### Phenotype Data Integration Framework

```mermaid
graph TD
    A[Multiple Phenotype Sources] --> B[Data Harmonization]
    B --> C[Standard Schema]

    A --> D[Source-specific Processing]
    D --> E[Format Conversion]
    E --> F[Unit Standardization]

    C --> G[Integrated Dataset]
    F --> G

    G --> H[Quality Assessment]
    H --> I[Completeness Check]
    I --> J[Consistency Check]
    J --> K[Outlier Detection]

    K --> L{Clean Data?}
    L -->|Yes| M[Data Cleaning]
    L -->|No| N[Flag Issues]

    M --> O[Final Dataset]
    N --> O

    O --> P[Cross-source Validation]
    P --> Q[Integrated Phenotypes]

    Q --> R[Downstream Analysis]
    R --> S[Biological Insights]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style G fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style S fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Data Sources"
        T[AntWiki] -.-> A
        U[Field Observations] -.-> A
        V[Literature] -.-> A
        W[Databases] -.-> A
    end

    subgraph "Harmonization Steps"
        X[Term Mapping] -.-> B
        Y[Unit Conversion] -.-> F
        Z[Scale Normalization] -.-> F
    end

    subgraph "Quality Metrics"
        AA[Missing Rate] -.-> H
        BB[Measurement Error] -.-> H
        CC[Biological Plausibility] -.-> H
    end
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

### AntWiki Web Scraping (`scraper.py`)
Comprehensive web scraping of AntWiki species pages to extract all sections including measurements, traits, taxonomy, distribution, and descriptions.

**Cloudflare Protection:** The scraper automatically uses `cloudscraper` if available to bypass Cloudflare protection. Install with:
```bash
pip install cloudscraper
# Or install optional dependencies:
pip install metainformant[scraping]
```

If `cloudscraper` is not available, the scraper falls back to `requests` but may fail on Cloudflare-protected sites.

**Usage:**
```python
from metainformant.phenotype.scraper import AntWikiScraper, load_scraper_config
from pathlib import Path

# Load configuration (with environment variable overrides)
config = load_scraper_config()

# Initialize scraper (automatically uses cloudscraper if available)
scraper = AntWikiScraper(config)

# Scrape single species
data = scraper.scrape_species_page("Camponotus_pennsylvanicus")
print(f"Found {len(data['traits'])} traits, {len(data['measurements'])} measurements")

# Scrape all species (with limit for testing)
stats = scraper.scrape_all_species(output_dir=Path("output/phenotype/antwiki/"), limit=10)
print(f"Completed: {stats['completed']}, Failed: {stats['failed']}")
```

**Command-line Usage:**
```bash
# Scrape single species
python3 scripts/phenotype/scrape_antwiki.py --species Camponotus_pennsylvanicus

# Scrape all species with limit
python3 scripts/phenotype/scrape_antwiki.py --all --limit 10 --delay 2.0

# Resume interrupted scraping
python3 scripts/phenotype/scrape_antwiki.py --all --resume
```

**Configuration:**
Scraper configuration is loaded from `config/phenotype/antwiki_scraper.yaml` with environment variable overrides:
- `PHEN_SCRAPE_DELAY`: Override delay between requests (seconds)
- `PHEN_WORK_DIR`: Override output directory

**Features:**
- **Cloudflare bypass**: Automatically uses `cloudscraper` if available to handle Cloudflare protection
- Rate limiting with randomized delays (50-150% of base delay) to mimic human behavior
- Robots.txt compliance checking
- Retry logic with exponential backoff
- Progress tracking and checkpoint/resume support
- Comprehensive data extraction from all page sections
- Browser-like headers including Referer for legitimate-looking requests
- Organized output structure (individual species files + combined dataset)

**Output Structure:**
```
output/phenotype/antwiki/
├── species/
│   ├── Camponotus_pennsylvanicus.json
│   └── ...
├── all_species.json
├── scraping_log.jsonl
└── checkpoint.json
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

## See Also

- **[AGENTS.md](AGENTS.md)**: AI agent contributions and development details for the phenotype module

## Related Modules

The Phenotype module integrates with several other METAINFORMANT modules:

- **Life Events Module**: Temporal phenotype analysis and life course event integration; longitudinal trait studies
- **GWAS Module**: Phenotype data for genome-wide association studies; quantitative and qualitative trait analysis
- **Networks Module**: Phenotype correlation networks, trait co-occurrence analysis, and morphological relationship modeling
- **Multi-omics Module**: Phenotype integration with DNA, RNA, and protein data; systems-level trait analysis
- **ML Module**: Machine learning prediction of phenotypes from genomic data; trait classification models
- **Visualization Module**: Phenotype distribution plots, trait correlation heatmaps, and morphological data visualization
- **Quality Module**: Phenotype data validation and quality control; trait measurement standardization
- **Ontology Module**: Phenotype annotation using standardized vocabularies; trait classification systems
- **Information Module**: Information-theoretic analysis of phenotype complexity and trait relationships
- **Simulation Module**: Synthetic phenotype generation for statistical power analysis and method validation

This module provides essential tools for phenotype-genotype association studies and morphological analysis.
