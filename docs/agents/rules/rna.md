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
  deconvolution/  bulk RNA-seq cell-type deconvolution
  engine/         streaming producer, workflow, progress, provenance
  integration/    adapters where RNA code reaches other domains
  retrieval/      ENA/SRA retrieval helpers
  splicing/       alternative-splicing analysis
  steps.py        registry of workflow step runners (STEP_RUNNERS)
```

## Current execution contract

```text
metadata → select → getfastq → integrate → quant
         → merge → wsfilter → finalize → sanity
```

`cstmm` and `csfilter` are inserted after `wsfilter` only when ortholog inputs
(`orthogroup_table` or `dir_busco`) are configured.

The producer is `StreamingPipelineOrchestrator`; progress is stored in
`ProgressDB`; downstream checkpoints (`merge`, `wsfilter`, `finalize`,
`sanity`) are owned by the Hymenoptera project scripts after the producer lock
is released.

## Current commands

```bash
export AMALGKIT_DATA_ROOT=/Volumes/external_drive/Data/amalgkit
uv run python scripts/rna/validate_configs.py
uv run python scripts/rna/validate_all_species_workflow.py
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

`projects/hymenoptera_amalgkit` is a git submodule; the project YAMLs and
`run_full_campaign.sh` exist only once it is checked out.

For one species, use `scripts/rna/process_species.py --species <name>` with the
same config directory and data root. Never start it while the full campaign
lock (`$AMALGKIT_DATA_ROOT/results/.full_campaign.lock`) is held.

## Public Python interfaces

```python
from pathlib import Path

from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator
from metainformant.rna.engine.workflow import load_workflow_config, plan_workflow

data_root = Path("/Volumes/external_drive/Data/amalgkit")
config = load_workflow_config(
    "projects/hymenoptera_amalgkit/config/amalgkit/amalgkit_apis_mellifera.yaml"
)
assert [name for name, _ in plan_workflow(config)][0] == "metadata"
db = ProgressDB(data_root / "pipeline_progress.db")
```

## Configuration

Use the canonical project YAMLs and shared selection rule file
(`select_rules.tsv`). Runtime budgets use the current `AMALGKIT_PIPELINE_*`
namespace (for example `AMALGKIT_PIPELINE_QUANT_SLOTS` and
`AMALGKIT_PIPELINE_MAX_IN_FLIGHT`) and the CLI options shown by `--help`;
`--max-gb`, `--workers` and `--threads` take their defaults from
`PIPELINE_MAX_GB`, `PIPELINE_WORKERS` and `PIPELINE_THREADS` when those are set.
Keep large artifacts under `$AMALGKIT_DATA_ROOT`; without it the data root is
`output/amalgkit` relative to the current directory.

## Idempotence and evidence

- Discovery and dry-run commands are read-only and deterministic.
- SQLite sample state is the source for resume and status.
- Reuse requires a hash-bound quantification provenance sidecar whose contract
  verifies: the accession matches and the SHA-256 of the configuration and
  reference manifest match the recorded values. The quantification file's own
  SHA-256 is rechecked on resume only when `AMALGKIT_RESUME_VERIFY_HASHES` is
  set. Outputs without a verifiable sidecar are quarantined.
  `--requantification-policy` chooses which verified outputs are re-quantified:
  `preserve` keeps all of them, `version-drift` re-runs those recorded under a
  different Amalgkit runtime, and `all` re-runs every one.
- Downstream stages validate exact inputs and write atomic outputs.
- Reports preserve failed, missing, and incomplete species.
- Cross-species analysis admits only finalized species that pass matrix and
  provenance gates.

## Testing

Use `uv run pytest` with real implementations. Mark network and external-tool
tests explicitly (`@pytest.mark.network`, `@pytest.mark.external_tool`), keep
test artifacts under `tmp_path`, run the project docs validator
(`projects/hymenoptera_amalgkit/scripts/validate_project_docs.py`), and run the
RNA current-surface retirement guard (`tests/rna/test_rna_amalgkit_current.py`)
before release. The guard fails when a retired entry-point file reappears or an
active RNA source, script, doc or config file (this page included) refers to a
retired API.
