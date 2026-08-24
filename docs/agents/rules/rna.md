# Current RNA module rules

## Purpose

RNA transcriptomic acquisition, quantification, evidence validation, and
analysis through the pinned Amalgkit 0.16.60 contract.

## Source structure

```text
src/metainformant/rna/
  amalgkit/       command registry and wrappers
  analysis/       expression, QC, cross-species, validation
  core/           configuration, cleanup, dependencies, environment
  engine/         streaming producer, workflow, progress, provenance
  retrieval/      ENA/SRA retrieval helpers
  splicing/       alternative-splicing analysis
```

## Current execution contract

```text
metadata → select → getfastq → integrate → quant
         → merge → wsfilter → finalize → sanity
```

The producer is `StreamingPipelineOrchestrator`; progress is stored in
`ProgressDB`; downstream checkpoints are owned by the Hymenoptera project
scripts after the producer lock is released.

## Current commands

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
uv run python scripts/rna/validate_configs.py
uv run python scripts/rna/validate_all_species_workflow.py
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

For one species, use `scripts/rna/process_species.py` with the same config
directory and data root. Never start it while the full campaign lock is held.

## Public Python interfaces

```python
from pathlib import Path

from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator
from metainformant.rna.engine.workflow import load_workflow_config, plan_workflow

data_root = Path("/Volumes/blue/data/amalgkit")
config = load_workflow_config(
    "projects/hymenoptera_amalgkit/config/amalgkit/amalgkit_apis_mellifera.yaml"
)
assert [name for name, _ in plan_workflow(config)][0] == "metadata"
db = ProgressDB(data_root / "pipeline_progress.db")
```

## Configuration

Use the canonical project YAMLs and shared selection rule file. Runtime
budgets use the current `AMALGKIT_PIPELINE_*` namespace and the CLI options
shown by `--help`. Keep large artifacts under `$AMALGKIT_DATA_ROOT`.

## Idempotence and evidence

- Discovery and dry-run commands are read-only and deterministic.
- SQLite sample state is the source for resume and status.
- Reuse requires current, hash-bound quantification provenance.
- Downstream stages validate exact inputs and write atomic outputs.
- Reports preserve failed, missing, and incomplete species.
- Cross-species analysis admits only finalized species that pass matrix and
  provenance gates.

## Testing

Use `uv run pytest` with real implementations. Mark network and external-tool
tests explicitly, keep test artifacts under `tmp_path`, run the project docs
validator, and run the RNA current-surface retirement guard before release.
