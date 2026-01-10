# Amalgkit

## Overview
This module provides a **robust, production-grade wrapper** for `amalgkit` within the MetaInformAnt framework. It orchestrates the full RNA-seq lifecycle: metadata retrieval, robust SRA downloading, automated genome indexing, quantification, and downstream analysis.

## Core Capabilities (Updated Jan 2026)
*   **Robust SRA Download**: Parallel, resumable downloads via AWS Open Data, bypassing fragile `prefetch` logic for LITE files.
*   **Automated Genome Prep**: Auto-fetches transcriptomes and builds kallisto indices if missing.
*   **Intelligent Workflow**: 
    *   **Smart Skipping**: Defer to per-sample incomplete status rather than skipping entire phases.
    *   **Safety Guards**: "Intelligent Redo" prevents accidental deletion of valid SRA files.
*   **Type-Safe Config**: Dataclass-based configuration with rigorous validation.
*   **No-Mock Testing**: All workflows verified against real external tools (SRA Toolkit, Kallisto, Amalgkit).

## Usage
Import and use the orchestration functions from `metainformant.rna.engine.workflow`:

```python
from metainformant.rna.engine.workflow import load_workflow_config, execute_workflow

config = load_workflow_config("config/my_species.yaml")
execute_workflow(config, steps=['getfastq', 'quant'])
```
