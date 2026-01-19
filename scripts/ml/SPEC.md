# SPEC: ML Scripts

Machine learning pipelines for biological trait prediction and feature selection.

## Workflows

- `train_trait_classifier.py`: Trains and cross-validates genomic or phenotypic classifiers.
- `select_important_features.py`: Uses mutual information or PCA to Identify key drivers of variation.

## Standards

- **Reproducibility**: Always log random seeds and hyperparameter configurations.
- **Evaluation**: Report ROC-AUC, F1-score, and Confusion Matrices for all classification tasks.
