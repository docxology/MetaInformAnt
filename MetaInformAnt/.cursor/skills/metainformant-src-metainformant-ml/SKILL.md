---
name: metainformant-src-metainformant-ml
description: METAINFORMANT rules for directory src/metainformant/ml. Use when editing, adding tests, or reviewing code under this path. Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations.
---

# METAINFORMANT — `src/metainformant/ml`

Before editing files in this subtree:

- Read [`AGENTS.md`](../../../src/metainformant/ml/AGENTS.md) for this folder (canonical technical context).
- Optional overview: [`README.md`](../../../src/metainformant/ml/README.md).
- Global rules: [`CLAUDE.md`](../../../CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, real implementations).
- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`](../../../docs/REAL_IMPLEMENTATION_POLICY.md).
- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.

## Module surface (generated, validated)
- Purpose: Machine learning module for METAINFORMANT.

This module provides machine learning capabilities for bioinformatics analysis,
including classification, regression, dimensionality reduction, feature
selection/engineering, and local LLM inference.
- Public submodules: `automl`, `deep_learning`, `evaluation`, `features`, `interpretability`, `llm`, `models`.
- Canonical import: `import metainformant.ml` (submodules: `from metainformant import ml` then `ml.<submodule>`).
- Test entry point: `uv run pytest tests/ml -q` (one pytest directory per invocation).

Keep changes scoped; match existing patterns in this directory.
