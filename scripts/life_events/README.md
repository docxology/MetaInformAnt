# Life Events Analysis Scripts

Life course event sequence analysis, embedding learning, and outcome prediction workflow orchestrators.

## Directory Structure

```
scripts/life_events/
├── run_life_events_analysis.py     # Main life events analysis orchestrator ⭐
├── generate_synthetic_data.py      # Synthetic life event data generation
├── learn_embeddings.py             # Sequence embedding learning
├── train_model.py                  # Outcome prediction model training
├── predict_outcomes.py             # Outcome prediction and interpretation
├── generate_visualizations.py  # Visualization suite
├── generate_statistical_summary.py # Statistical analysis and reporting
├── life_course_example.py          # Life course analysis examples
├── analyze_intervention.py         # Intervention impact analysis
├── compare_groups.py               # Group comparison utilities
├── export_embeddings.py            # Embedding export and serialization
├── interpret_predictions.py        # Prediction interpretation tools
├── validate_data.py                # Data validation and quality checks
├── visualize_sequences.py          # Sequence visualization utilities
└── README.md                       # This file
```

## Main Analysis Workflow (`run_life_events_analysis.py`)

Comprehensive life course analysis orchestrator combining data generation, embedding learning, model training, and visualization.

**Features:**
- Synthetic life event sequence generation
- Sequence embedding learning (Word2Vec/Skip-gram)
- Outcome prediction modeling
- Intervention analysis
- Comprehensive statistical reporting
- Complete visualization suite

**Usage:**
```bash
# Generate synthetic data and run complete analysis
python3 scripts/life_events/run_life_events_analysis.py --synthetic --n-sequences 100 --generate-outcomes

# Analyze existing sequences with config file
python3 scripts/life_events/run_life_events_analysis.py --input data/life_events/sequences.json --config config/life_events_template.yaml

# Run with custom embedding parameters
python3 scripts/life_events/run_life_events_analysis.py --synthetic --embedding-dim 200 --window-size 10 --epochs 20

# Generate visualizations only
python3 scripts/life_events/run_life_events_analysis.py --input data/life_events/sequences.json --visualize-only
```

## Data Generation and Processing

### Synthetic Data Generation (`generate_synthetic_data.py`)

Generate realistic synthetic life event sequences for testing and validation.

**Features:**
- Configurable sequence length and complexity
- Realistic event transitions and probabilities
- Demographic and outcome variable generation
- Temporal dependency modeling

**Usage:**
```bash
# Generate synthetic life event data
python3 scripts/life_events/generate_synthetic_data.py --n-sequences 1000 --max-length 50 --output output/life_events/synthetic
```

### Data Validation (`validate_data.py`)

Validate life event sequence data quality and completeness.

**Usage:**
```bash
# Validate life event data
python3 scripts/life_events/validate_data.py --input data/life_events/sequences.json
```

## Machine Learning and Analysis

### Embedding Learning (`learn_embeddings.py`)

Learn vector representations of life event sequences using neural language models.

**Features:**
- Word2Vec/Skip-gram architectures
- Configurable embedding dimensions
- Sequence context window tuning
- Model serialization and export

**Usage:**
```bash
# Learn sequence embeddings
python3 scripts/life_events/learn_embeddings.py --input data/life_events/sequences.json --embedding-dim 150 --output output/life_events/embeddings
```

### Model Training (`train_model.py`)

Train predictive models for life course outcomes using learned embeddings.

**Features:**
- Multiple classifier architectures
- Cross-validation and hyperparameter tuning
- Feature importance analysis
- Model performance metrics

**Usage:**
```bash
# Train outcome prediction model
python3 scripts/life_events/train_model.py --embeddings output/life_events/embeddings.npy --outcomes data/life_events/outcomes.csv --output output/life_events/models
```

### Outcome Prediction (`predict_outcomes.py`)

Generate predictions and probability estimates for life course outcomes.

**Usage:**
```bash
# Generate outcome predictions
python3 scripts/life_events/predict_outcomes.py --model output/life_events/models --sequences data/life_events/test_sequences.json --output output/life_events/predictions
```

## Analysis and Interpretation

### Statistical Summary (`generate_statistical_summary.py`)

Comprehensive statistical analysis of life event sequences and outcomes.

**Features:**
- Sequence length and complexity statistics
- Event frequency and transition analysis
- Outcome distribution analysis
- Temporal pattern identification

**Usage:**
```bash
# Generate statistical summary
python3 scripts/life_events/generate_statistical_summary.py --input data/life_events/sequences.json --output output/life_events/statistics
```

### Intervention Analysis (`analyze_intervention.py`)

Analyze the impact of interventions on life course trajectories.

**Usage:**
```bash
# Analyze intervention effects
python3 scripts/life_events/analyze_intervention.py --sequences data/life_events/sequences.json --interventions data/interventions.json --output output/life_events/intervention_analysis
```

### Group Comparison (`compare_groups.py`)

Compare life event patterns between different demographic or outcome groups.

**Usage:**
```bash
# Compare groups
python3 scripts/life_events/compare_groups.py --sequences data/life_events/sequences.json --groups data/groups.csv --output output/life_events/group_comparison
```

## Visualization

### Visualization Suite (`generate_visualizations.py`)

Generate comprehensive visualizations for life event analysis results.

**Features:**
- Sequence pattern visualizations
- Embedding space projections
- Model performance plots
- Outcome prediction visualizations
- Statistical analysis plots

**Usage:**
```bash
# Generate all visualizations
python3 scripts/life_events/generate_visualizations.py --results-dir output/life_events --output output/life_events/visualizations
```

### Sequence Visualization (`visualize_sequences.py`)

Visualize individual and aggregate life event sequences.

**Usage:**
```bash
# Visualize sequences
python3 scripts/life_events/visualize_sequences.py --sequences data/life_events/sequences.json --output output/life_events/plots
```

## Utilities

### Embedding Export (`export_embeddings.py`)

Export learned embeddings in various formats for downstream analysis.

**Usage:**
```bash
# Export embeddings
python3 scripts/life_events/export_embeddings.py --model output/life_events/embeddings.model --format tsv --output output/life_events/embeddings_export
```

### Prediction Interpretation (`interpret_predictions.py`)

Interpret and explain model predictions with feature importance and decision paths.

**Usage:**
```bash
# Interpret predictions
python3 scripts/life_events/interpret_predictions.py --predictions output/life_events/predictions.json --model output/life_events/models --output output/life_events/interpretation
```

### Life Course Examples (`life_course_example.py`)

Example workflows demonstrating life course analysis techniques.

**Usage:**
```bash
# Run life course examples
python3 scripts/life_events/life_course_example.py --example basic-sequence-analysis
```

## Output Structure

```
output/life_events/
├── synthetic_data/                # Generated synthetic sequences
│   ├── sequences.json
│   ├── outcomes.csv
│   └── metadata.json
├── embeddings/                    # Learned sequence embeddings
│   ├── model.pkl
│   ├── embeddings.npy
│   └── vocabulary.json
├── models/                        # Trained predictive models
│   ├── classifier.pkl
│   ├── feature_importance.json
│   └── performance_metrics.json
├── predictions/                   # Outcome predictions
│   ├── probabilities.csv
│   ├── classifications.csv
│   └── confidence_scores.json
├── statistics/                    # Statistical summaries
│   ├── sequence_stats.json
│   ├── transition_matrices.json
│   └── outcome_distributions.json
├── visualizations/                # All generated plots
│   ├── sequence_patterns.png
│   ├── embedding_projections.png
│   ├── model_performance.png
│   ├── outcome_predictions.png
│   └── statistical_analysis.png
└── analysis_report.json           # Comprehensive analysis report
```

## Key Features

✅ **End-to-End Workflow**: Data generation → embedding → training → prediction → visualization
✅ **Synthetic Data Generation**: Realistic life event sequences for testing
✅ **Advanced ML**: Word2Vec embeddings and predictive modeling
✅ **Comprehensive Analysis**: Statistical summaries and intervention analysis
✅ **Rich Visualization**: Complete plotting suite for all analysis types
✅ **Flexible Input**: Support for both synthetic and real life event data
✅ **Modular Design**: Individual scripts for specific analysis components

## Integration

Integrates with:
- **metainformant.life_events**: Core life course analysis functionality
- **Machine learning**: scikit-learn, gensim for embeddings
- **Statistical analysis**: pandas, numpy, scipy
- **Visualization**: matplotlib, seaborn, plotly
- **Core utilities**: I/O, logging, configuration management

## Dependencies

- **metainformant.life_events**: Life course analysis module
- **gensim**: Word2Vec embeddings
- **scikit-learn**: Machine learning algorithms
- **pandas/numpy**: Data manipulation and statistics
- **matplotlib/seaborn**: Visualization support
- **plotly**: Interactive visualizations

## Related Documentation

- [Life Events Analysis Documentation](../../docs/life_events/README.md)
- [Machine Learning Integration](../../docs/life_events/ml_integration.md)
- [Synthetic Data Generation](../../docs/life_events/synthetic_data.md)
- [METAINFORMANT CLI](../../docs/cli.md)

