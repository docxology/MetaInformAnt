# Current RNA API

The API below is the supported Python surface for the current Amalgkit
workflow. The command registry and project scripts remain the authoritative
execution contracts.

## Streaming producer

```python
from pathlib import Path

from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator

data_root = Path("/Volumes/blue/data/amalgkit")
producer = StreamingPipelineOrchestrator(
    config_dir=Path("projects/hymenoptera_amalgkit/config/amalgkit"),
    log_dir=data_root / "logs",
    db_path=data_root / "pipeline_progress.db",
)
producer.run_all(
    ["amalgkit_apis_mellifera.yaml"],
    max_gb=50.0,
    workers=4,
    threads=8,
    quant_slots=2,
)
```

`run_all()` handles metadata, selection, acquisition, integration, and
quantification. It uses `ProgressDB` for resumable state and writes current
provenance for every reusable quantification.

## Progress database

```python
from metainformant.rna.engine.progress_db import ProgressDB

db = ProgressDB(data_root / "pipeline_progress.db")
print(db.get_counts())
```

The sample states are `pending`, `downloading`, `downloaded`, `quantifying`,
`quantified`, `quarantined`, and `failed`. The database is the source for
progress reporting; directory scans are not completion evidence. Quantification
audit records distinguish current, compatible version drift, legacy-unverified,
and invalid outputs.

## Workflow planning

```python
from metainformant.rna.engine.workflow import load_workflow_config, plan_workflow

config = load_workflow_config(
    "projects/hymenoptera_amalgkit/config/amalgkit/amalgkit_apis_mellifera.yaml"
)
steps = plan_workflow(config)
assert [name for name, _ in steps] == [
    "metadata", "select", "getfastq", "integrate", "quant",
    "merge", "wsfilter", "finalize", "sanity",
]
```

(quant)=
## Quantification API

The quantification wrapper runs the pinned Amalgkit quantification command and
returns a `subprocess.CompletedProcess`. A successful process is reusable only
when the accompanying hash-bound quantification provenance receipt is current.

```python
from metainformant.rna.amalgkit import quant

result = quant(
    out_dir="output/amalgkit/work",
    metadata="output/amalgkit/work/metadata/pivot_qualified.tsv",
    index_dir="output/amalgkit/work/index",
    threads=8,
    clean_fastq=False,
)
```

## Provenance

The `metainformant.rna.engine.provenance` module writes and validates current
hash-bound receipts for metadata, references, quantification, downstream
stages, and raw-input reclamation. A missing or mismatched receipt causes the
corresponding output to be treated as incomplete.

## Analysis namespaces

`metainformant.rna.analysis.expression`, `.qc`, `.cross_species`, and
`.validation` provide table-level analysis helpers. They consume finalized,
evidence-checked inputs; they do not promote incomplete production outputs to
biological results.
