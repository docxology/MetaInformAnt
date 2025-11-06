# Life Events Module

The life events module provides comprehensive tools for analyzing human life courses as temporal event sequences, enabling prediction of life outcomes from event data using NLP-inspired methods.

## Overview

This module treats life courses as sequences of events (similar to how NLP treats text as sequences of words), allowing the application of sequence modeling techniques to predict outcomes like early mortality, personality traits, and health outcomes.

## Key Features

- **Event Sequence Representation**: Structured data types for temporal event sequences with domain categorization
- **Event Embeddings**: Learn dense vector representations of events using Word2Vec-style methods
- **Sequence Embeddings**: Aggregate event sequences into fixed-size vectors for ML models
- **Domain-Specific Embeddings**: Separate embedding spaces for different life domains
- **Prediction Models**: Classification and regression models for predicting life outcomes
- **Visualization**: Timeline plots, embedding visualizations, attention heatmaps
- **Model Interpretation**: Event importance, temporal patterns, feature attribution
- **Workflow Integration**: End-to-end pipelines for complete analysis

## Documentation

- **[Complete Module Documentation](README.md)** - Comprehensive guide with examples and API reference
- **[AI Contributions](AGENTS.md)** - AI assistance in module development

## Quick Start

```python
from metainformant.life_events import Event, EventSequence, analyze_life_course
from datetime import datetime

# Create event sequences
events = [
    Event("degree", datetime(2010, 6, 1), "education"),
    Event("job_change", datetime(2015, 3, 1), "occupation"),
]
sequence = EventSequence(person_id="person_001", events=events)

# Analyze life course
results = analyze_life_course([sequence], outcomes=None)
```

## Core Components

### Event Data Structures
- `Event`: Individual life event with temporal and domain information
- `EventSequence`: Container for temporal event sequence of a single person
- `EventDatabase`: Collection of multiple event sequences

### Event Embeddings
- Word2Vec-style embedding learning (Skip-gram, CBOW)
- Domain-specific embedding spaces
- Sequence embedding aggregation methods

### Prediction Models
- Classification and regression models
- LSTM and GRU sequence models
- Ensemble and multi-task predictors
- Survival analysis models

### Visualization
- Timeline plots
- Embedding visualizations
- Attention heatmaps
- Importance plots

## Integration

The life events module integrates with:
- **ML**: Classification and regression infrastructure
- **Phenotype**: Life course phenotype extraction
- **Visualization**: Plotting infrastructure
- **Multi-Omics**: Event sequences as one data layer

## Related Documentation

- **[Source Module Documentation](../src/metainformant/life_events/README.md)** - Detailed implementation
- **[ML Module](../ml/index.md)** - Machine learning methods
- **[Phenotype Module](../phenotype/index.md)** - Phenotypic trait analysis
- **[Visualization Module](../visualization/index.md)** - Plotting utilities

