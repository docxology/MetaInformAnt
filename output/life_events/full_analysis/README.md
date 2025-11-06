# Life Events Full Analysis Workflow Results

## Overview

This directory contains the complete results from running the full life events analysis workflow with synthetic data generation, statistical analysis, and comprehensive visualization.

## Generated Data

### Input Data
- **Synthetic Sequences**: `synthetic_sequences.json` (237 KB)
  - 100 event sequences
  - 1,752 total events
  - 42 unique event types
  - Mean sequence length: 17.5 events
  - Temporal span: ~8.9 years (mean)

- **Synthetic Outcomes**: `synthetic_outcomes.json` (534 B)
  - 100 outcome values
  - Generated with complex relationship to events

### Analysis Results
- **Event Embeddings**: `embeddings.json` (91 KB)
  - Embeddings for 42 unique events
  - 100-dimensional embeddings
  - Learned using Skip-gram with window size 5

- **Prediction Model**: `model.json` (414 KB)
  - Embedding-based regression model
  - Trained on 100 sequences
  - 15 training epochs

- **Analysis Results**: `analysis_results.json` (2.3 KB)
  - Model predictions for all sequences
  - Performance metrics

- **Statistical Summary**: `statistical_summary.json` (3.7 KB)
  - Comprehensive sequence statistics
  - Temporal pattern analysis
  - Event co-occurrence patterns
  - Outcome statistics
  - Model performance metrics

- **Workflow Summary**: `workflow_summary.json` (3.1 KB)
  - Complete workflow configuration
  - Final results summary

## Visualizations (23 files)

### Sequence Analysis
1. **Event Timeline** (`event_timeline_first_sequence.png`)
   - Visual timeline of events in first sequence
2. **Timeline Sequences** (`timeline_sequence_*.png`)
   - Individual timeline visualizations for multiple sequences
3. **Sequence Length Distribution** (`sequence_length_distribution.png`)
   - Histogram of sequence lengths

### Domain Analysis
4. **Domain Distribution** (`domain_distribution.png`)
   - Bar/pie chart of event frequencies by domain
5. **Domain Timeline** (`domain_timeline.png`)
   - Gantt-style timeline across multiple domains

### Temporal Analysis
6. **Temporal Density** (`temporal_density.png`)
   - Histogram of event density over time
7. **Temporal Patterns** (`temporal_patterns.png`)
   - Time-based importance visualization
8. **Event Frequency Heatmap** (`event_frequency_heatmap.png`)
   - Temporal frequency heatmap by domain

### Event Relationships
9. **Event Co-occurrence** (`event_cooccurrence.png`)
   - Heatmap of event co-occurrence patterns
10. **Transition Network** (`transition_network.png`)
    - Network graph of event transitions
11. **Sequence Similarity** (`sequence_similarity.png`)
    - Cosine similarity matrix between sequences

### Embeddings
12. **Event Embeddings PCA** (`event_embeddings_pca.png`)
    - 2D PCA projection of event embeddings
13. **Embedding Clusters** (`embedding_clusters.png`)
    - Clustered visualization of embeddings (UMAP)

### Outcomes & Predictions
14. **Outcome Distribution** (`outcome_distribution.png`)
    - Histogram/box plot of outcome values
15. **Prediction Accuracy** (`prediction_accuracy.png`)
    - Confusion matrix/ROC curve for predictions
16. **Prediction Importance** (`prediction_importance.png`)
    - Feature importance for predictions

## Key Statistics

### Sequence Statistics
- **Total Sequences**: 100
- **Total Events**: 1,752
- **Unique Event Types**: 42
- **Mean Sequence Length**: 17.5 events
- **Mean Temporal Span**: 3,236 days (~8.9 years)

### Domain Distribution
- **Health**: 703 events (40.1%)
- **Education**: 375 events (21.4%)
- **Occupation**: 233 events (13.3%)
- **Income**: 180 events (10.3%)
- **Address**: 124 events (7.1%)
- **Other**: 137 events (7.8%)

### Temporal Patterns
- **Date Range**: 2015-11-09 to 2025-11-03 (10 years)
- **Mean Inter-Event Interval**: 195.4 days
- **Median Inter-Event Interval**: 132 days
- **Events by Day of Week**: Relatively uniform distribution
- **Events by Month**: Relatively uniform distribution

### Top Event Types
1. Health: recovery (98 events)
2. Health: hospitalization (92 events)
3. Health: surgery (93 events)
4. Health: checkup (91 events)
5. Health: vaccination (88 events)
6. Education: graduation (61 events)
7. Education: degree (57 events)
8. Education: training (57 events)
9. Education: certification (55 events)

### Event Co-occurrences
Top co-occurring event pairs:
- Health: recovery | Health: vaccination (96 co-occurrences)
- Health: hospitalization | Health: surgery (104 co-occurrences)
- Health: checkup | Health: therapy (87 co-occurrences)
- Health: checkup | Health: vaccination (88 co-occurrences)

### Model Performance
- **Model Type**: Embedding-based regression
- **Embedding Dimension**: 100
- **Training Epochs**: 15
- **RMSE**: 0.00013
- **MAE**: 0.00010
- **Predictions Mean**: 1.000
- **Predictions Std**: 0.00013

## Workflow Execution

The complete workflow was executed using:

```bash
python3 scripts/life_events/run_life_events_analysis.py \
  --synthetic \
  --n-sequences 100 \
  --min-events 10 \
  --max-events 25 \
  --generate-outcomes \
  --outcome-relationship complex \
  --embedding-dim 100 \
  --window-size 5 \
  --epochs 15 \
  --model-type embedding \
  --all-visualizations \
  --output output/life_events/full_analysis \
  --verbose
```

Additional visualizations were generated using:
- `scripts/life_events/visualize_sequences.py` - Sequence visualizations
- `scripts/life_events/generate_all_visualizations.py` - Comprehensive visualization suite
- `scripts/life_events/generate_statistical_summary.py` - Statistical analysis

## Files Generated

**Total**: 30 files
- 7 JSON data/result files
- 23 PNG visualization files

**Total Size**: 8.8 MB

## Next Steps

1. **Explore Visualizations**: Review all PNG files in `visualizations/` directory
2. **Analyze Statistics**: Examine `statistical_summary.json` for detailed patterns
3. **Model Evaluation**: Review `analysis_results.json` for prediction performance
4. **Custom Analysis**: Use modular scripts in `scripts/life_events/` for specific analyses

## Scripts Used

- `run_life_events_analysis.py` - Main orchestrator
- `visualize_sequences.py` - Sequence visualizations
- `generate_all_visualizations.py` - Comprehensive visualization suite
- `generate_statistical_summary.py` - Statistical analysis

All scripts are located in `scripts/life_events/`.

---

*Generated: 2025-11-05*
*Analysis Type: Full Life Events Workflow*
*Dataset: Synthetic (100 sequences, 1,752 events)*

