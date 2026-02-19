# Agent Directives: ml

**Context**: Machine learning module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `automl/` — exports: `optimization`
- `deep_learning/`
- `evaluation/` — exports: `validation`
- `features/` — exports: `dimensionality`, `features`
- `interpretability/` — exports: `explainers`, `feature_selection`
- `llm/` — exports: `ollama`
- `models/` — exports: `classification`, `regression`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
