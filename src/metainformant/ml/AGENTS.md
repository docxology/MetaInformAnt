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
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow REAL IMPLEMENTATION policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/ml/](../../../docs/ml/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Core module**: [../core/AGENTS.md](../core/AGENTS.md)
- **Gwas module**: [../gwas/AGENTS.md](../gwas/AGENTS.md)
- **Networks module**: [../networks/AGENTS.md](../networks/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
