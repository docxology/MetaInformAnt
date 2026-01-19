# SPEC: Life Events Scripts

Scripts for processing temporal event sequences and training embedding models.

## Workflows

- `analyze_life_events.py`: Basic statistical analysis of event sequences.
- `train_embeddings.py`: Orchestrates the training of Word2Vec or LSTM models on life events data.

## Standards

- **Anonymization**: All person-level data must be anonymized before processing.
- **Reporting**: Output model performance metrics to `output/life_events/metrics/`.
