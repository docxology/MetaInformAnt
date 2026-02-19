# Agent Directives: life_events

**Context**: Life events and trajectory analysis module for METAINFORMANT.



## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `analysis/` — exports: `interpretability`
- `core/` — exports: `config`, `events`, `utils`
- `models/` — exports: `embeddings`, `predictor`, `sequence_models`, `statistical_models`
- `survival/` — exports: `time_to_event`
- `visualization/` — exports: `network`, `statistical`, `timeline`
- `workflow/` — exports: `workflow`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
