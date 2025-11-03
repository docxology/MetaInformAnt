# Ecology Module

The `ecology` module provides tools for ecological metadata management, community analysis, and environmental data integration.

## Overview

This module handles ecological data including species diversity, community composition, and environmental parameters that influence biological systems.

## Submodules

### Community Analysis (`community.py`)
Tools for analyzing ecological communities and species interactions.

**Key Features:**
- Species diversity metrics (Shannon, Simpson indices)
- Community composition analysis
- Species abundance distributions
- Beta diversity calculations

**Usage:**
```python
from metainformant.ecology.community import (
    shannon_diversity,
    simpson_diversity,
    species_richness,
    pielou_evenness,
    chao1_estimator
)

# Community analysis
abundance_data = [10, 8, 6, 4, 2, 1]
diversity = shannon_diversity(abundance_data)
simpson = simpson_diversity(abundance_data)
richness = species_richness(abundance_data)
evenness = pielou_evenness(abundance_data)
richness_est = chao1_estimator([int(x) for x in abundance_data])
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import population
from metainformant.ecology.community import shannon_diversity

# Ecological context for genetic diversity
genetic_diversity = population.nucleotide_diversity(sequences)
# Ecological diversity from species abundances
species_abundances = [10, 8, 6, 4, 2]
ecological_diversity = shannon_diversity(species_abundances)
```

### With Visualization Module
```python
from metainformant.ecology.community import community_metrics
from metainformant.visualization import heatmap
import numpy as np

# Visualize community structure
# Create abundance matrix (sites x species)
abundance_matrix = np.array([[10, 8, 6], [5, 7, 9], [12, 4, 8]])
ax = heatmap(abundance_matrix)
ax.set_title("Species Abundance")
```

## Data Sources

- Species occurrence databases
- Environmental monitoring data
- Biodiversity surveys and inventories
- Ecological metadata repositories

## Performance Features

- Efficient processing of large ecological datasets
- Memory-optimized community calculations
- Support for sparse ecological matrices

## Testing

Comprehensive tests cover:
- Diversity metric calculations
- Community analysis algorithms
- Integration with environmental data

## Dependencies

- NumPy for numerical computations
- Pandas for data manipulation
- Optional: specialized ecological packages

This module provides essential tools for ecological data analysis and biodiversity research.
