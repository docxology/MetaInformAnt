# Information Theory Module

The information theory module provides comprehensive information-theoretic methods for analyzing biological data, including both syntactic (Shannon entropy, mutual information) and semantic (information content, semantic similarity) information measures.

## Overview

This module implements fundamental information theory concepts for biological data analysis:
- **Syntactic Information Theory**: Shannon entropy, mutual information, KL divergence, transfer entropy
- **Semantic Information Theory**: Information content, semantic similarity, semantic entropy
- **Network Integration**: Information flow, network entropy, information-based community detection
- **Visualization**: Integration with visualization module for plotting information measures

## Documentation

- **[Complete Module Documentation](README.md)** - Comprehensive guide with examples and API reference
- **[AI Contributions](AGENTS.md)** - AI assistance in module development

## Quick Start

```python
from metainformant.information import shannon_entropy, mutual_information

# Calculate entropy of a probability distribution
probs = [0.5, 0.3, 0.2]
entropy = shannon_entropy(probs)  # bits

# Calculate mutual information between sequences
x = [0, 1, 0, 1, 0, 1]
y = [1, 0, 1, 0, 1, 0]
mi = mutual_information(x, y)
```

## Key Components

### Syntactic Information Theory
- Shannon entropy and variants (Rényi, Tsallis)
- Mutual information and normalized variants
- Kullback-Leibler and Jensen-Shannon divergence
- Transfer entropy for time series

### Semantic Information Theory
- Information content for hierarchical terms
- Semantic similarity measures
- Semantic entropy

### Network Integration
- Network entropy analysis
- Information flow through networks
- Information-based community detection

### Continuous Methods
- Differential entropy for continuous distributions
- Continuous mutual information estimation
- Bias correction methods

## Integration

The information module integrates with:
- **DNA**: Sequence information analysis
- **RNA**: Expression data entropy analysis
- **Networks**: Network information measures
- **Ontology**: GO term information content
- **Multi-Omics**: Cross-platform information analysis
- **ML**: Feature selection using mutual information

## Related Documentation

- **[Source Module Documentation](../src/metainformant/information/README.md)** - Detailed implementation
- **[Workflows Documentation](../src/metainformant/information/WORKFLOWS.md)** - Workflow guides
- **[Examples Documentation](../src/metainformant/information/EXAMPLES.md)** - Usage examples

