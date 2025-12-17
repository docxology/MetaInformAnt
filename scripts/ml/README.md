# Machine Learning Scripts

Machine learning pipeline orchestration for biological data analysis.

## Directory Structure

```
scripts/ml/
├── run_ml_pipeline.py             # Machine learning workflow orchestrator
└── README.md                      # This file
```

## Machine Learning Pipeline (`run_ml_pipeline.py`)

Comprehensive machine learning workflow orchestrator for classification, regression, feature selection, and model evaluation.

**Features:**
- Classification and regression tasks
- Feature selection and dimensionality reduction
- Cross-validation and hyperparameter tuning
- Model evaluation and performance metrics
- Feature importance analysis

**Usage:**
```bash
# Classification with feature selection
python3 scripts/ml/run_ml_pipeline.py --features X.csv --labels y.csv --classify --feature-selection --n-features 50

# Regression with cross-validation
python3 scripts/ml/run_ml_pipeline.py --features X.csv --labels y.csv --regress --cross-validate --cv-folds 5

# Dimensionality reduction
python3 scripts/ml/run_ml_pipeline.py --features X.csv --reduce-dimensions --method pca --n-components 20
```

**Options:**
- `--features`: Feature matrix file (CSV/TSV format)
- `--labels`: Target labels/values file
- `--classify`: Perform classification task
- `--regress`: Perform regression task
- `--feature-selection`: Enable feature selection
- `--cross-validate`: Perform cross-validation
- `--reduce-dimensions`: Apply dimensionality reduction
- `--output`: Output directory (defaults to output/ml/)
- `--threads`: Number of threads to use
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/ml/
├── feature_selection/             # Feature selection results
│   ├── selected_features.json
│   └── feature_importance.json
├── dimensionality_reduction/      # Dimension reduction results
│   ├── transformed_features.csv
│   └── explained_variance.json
├── model_training/                # Trained models
│   ├── classifier.pkl
│   ├── regressor.pkl
│   └── model_metadata.json
├── cross_validation/              # Cross-validation results
│   ├── cv_scores.json
│   └── performance_metrics.json
├── evaluation/                    # Model evaluation results
│   ├── classification_report.json
│   ├── confusion_matrix.json
│   └── roc_curves.json
└── analysis_report.json           # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.ml**: Core machine learning functionality
- **scikit-learn**: Machine learning algorithms and utilities
- **Core utilities**: I/O, logging, path management
- **Visualization**: Plot generation for results

## Dependencies

- **metainformant.ml**: Machine learning module
- **scikit-learn**: Machine learning algorithms
- **pandas/numpy**: Data manipulation
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Machine Learning Documentation](../../docs/ml/README.md)
- [Core Utilities](../../docs/core/README.md)
- [METAINFORMANT CLI](../../docs/cli.md)
