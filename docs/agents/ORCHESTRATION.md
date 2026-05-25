# Orchestration: Workflow Manager and Coordination Layer

**API Reference**: `metainformant.core.engine.workflow_manager`  
**Alternative**: `metainformant.core.execution.workflow` (config-driven)  
**Audience**: Developers building multi-agent pipelines, workflow engineers

## Table of Contents

- [Overview](#overview)
- [BasePipelineManager](#basepipelinemanager)
- [PipelinePhase](#pipelinephase)
- [PipelineItem](#pipelineitem)
- [Stage Lifecycle](#stage-lifecycle)
- [Phase Handler Patterns](#phase-handler-patterns)
- [Threading and Concurrency](#threading-and-concurrency)
- [TUI Integration](#tui-integration)
- [Configuration-Driven Workflows](#configuration-driven-workflows)
- [Advanced Techniques](#advanced-techniques)
- [Common Pitfalls](#common-pitfalls)

## Overview

The orchestration layer provides two primary APIs:

1. **Object-Oriented Pipeline** (`BasePipelineManager`) — for programmatic multi-phase workflows
2. **Declarative Configuration** (`download_and_process_data`) — for YAML/JSON-defined pipelines

Both converge on `metainformant.core.execution.parallel` for execution.

### When to Use Which?

| Scenario | Recommended API | Rationale |
|----------|-----------------|-----------|
| Fixed sequence of phases | `BasePipelineManager` | Explicit, type-safe, programmable logic |
| User-configurable pipeline | Config-driven (`workflow.py`) | End-users edit YAML; no code changes |
| Dynamically determined phases | `BasePipelineManager` (modify `self.phases` at runtime) | Conditional branching requires code |
| Simple download + process | Config-driven | Concise, declarative |
| Complex domain logic spanning 5+ phases | `BasePipelineManager` | Readable, testable handlers |

## BasePipelineManager

Generic, domain-agnostic pipeline orchestrator.

### Constructor

```python
from metainformant.core.engine.workflow_manager import BasePipelineManager, PipelinePhase

phases = [
    PipelinePhase("Download", download_handler, color=CYAN),
    PipelinePhase("Extract", extract_handler, color=MAGENTA),
    PipelinePhase("Quantify", quantify_handler, color=YELLOW),
]

manager = BasePipelineManager(
    phases=phases,
    max_threads=5,         # ThreadPoolExecutor size
    config=None,           # Optional dict of global config
)
```

**Parameters**:
- `phases`: Ordered list of `PipelinePhase` objects (executed left-to-right)
- `max_threads`: Worker pool size (default: 5). Use `resource_aware_workers()` for auto-tuning.
- `config`: Optional dict accessible to all handlers via `manager.config`

### Adding Work Items

```python
# Register a single item
item = manager.add_item(
    item_id="sample_001",
    metadata={"species": "Apis_mellifera", "accession": "SRR123"}
)

# Batch registration
for accession in sra_accessions:
    manager.add_item(accession, metadata={"source": "SRA"})
```

**Effect**: Adds item to `manager.items` dict, initializes TUI bar for that item.

### Running the Pipeline

```python
results = manager.run()
# Returns: Dict[str, bool] mapping item_id → success (True/False)
```

**Execution flow**:
1. Starts TUI (`manager.ui.start()`)
2. For each phase in `manager.phases`:
   - Filters eligible items via `phase.filter_fn`
   - Calls `phase.handler(manager, eligible_items)`
   - Handler must call `manager.mark_done()` or `manager.mark_failed()` for each item
3. Stops TUI, shuts down executor
4. Returns success dictionary

**Important**: If a handler raises an exception, pipeline aborts. Handlers must catch exceptions and call `mark_failed()` instead.

### State Inspection

During or after execution:

```python
# Access all items
for item_id, item in manager.items.items():
    print(f"{item_id}: {item.stage.value}")
    if item.error:
        print(f"  Error: {item.error}")

# Count by stage
from collections import Counter
stage_counts = Counter(item.stage for item in manager.items.values())
# e.g., Counter({'DONE': 150, 'FAILED': 5, 'PENDING': 0})
```

## PipelinePhase

Describes a single stage in the pipeline.

### Definition

```python
from metainformant.core.engine.workflow_manager import PipelinePhase

phase = PipelinePhase(
    name="Download",                    # Display name in TUI footer
    handler=download_handler,           # Callable: (manager, items) -> None
    filter_fn=lambda item: item.stage == Stage.PENDING,  # Eligibility predicate
    color=CYAN,                         # ANSI color for TUI
)
```

**Parameters**:
- `name`: Short label shown in UI footer and logs
- `handler`: Function that processes `items`. **Must** call `manager.mark_done()` or `manager.mark_failed()` for each item.
- `filter_fn`: Predicate selecting which items enter this phase. Default: `item.stage == Stage.PENDING`.
- `color`: TUI color constant (`BLUE`, `CYAN`, `GREEN`, `YELLOW`, `MAGENTA`, `RED`).

### Handler Signature

```python
def phase_handler(
    manager: BasePipelineManager,
    items: list[PipelineItem]
) -> None:
    """
    Process eligible items.

    Must call one of manager.mark_* for each item before returning.
    """
    for item in items:
        manager.mark_running(item)  # Optional but recommended
        # ... do work ...
        manager.mark_done(item)  # or manager.mark_failed(item, error_message)
```

**Concurrency**: Handler may submit work to `manager.executor` (ThreadPoolExecutor) for parallel execution. In that case, handler should return immediately after submission; completion callbacks call `mark_done()`/`mark_failed()`.

**Example**: See `_rna_download_phase()` in `workflow_manager.py` — submits parallel downloads, polls progress via file sizes, marks items upon `as_completed()`.

## PipelineItem

Tracked work item flowing through the pipeline.

### Structure

```python
from dataclasses import dataclass
from metainformant.core.engine.workflow_manager import PipelineItem, Stage

item = PipelineItem(
    item_id="SRR123",
    metadata={"species": "honeybee", "study": "Smith2023"},
    stage=Stage.PENDING,
    error=""
)
```

**Fields**:
- `item_id: str` — Unique identifier (e.g., SRA accession, sample barcode)
- `metadata: dict` — Arbitrary key-value bag passed between phases
- `stage: Stage` — Current lifecycle stage
- `error: str` — Human-readable error if `stage == FAILED`

### Metadata Propagation Pattern

Use `metadata` to pass computed values:

```python
# Phase 1: Download
item.metadata["sra_path"] = str(dest_path)
manager.mark_done(item)

# Phase 2: Extract (reads metadata from Phase 1)
sra_path = item.metadata["sra_path"]
fastq_dir = extract(sra_path)
item.metadata["fastq_dir"] = fastq_dir
manager.mark_done(item)

# Phase 3: Quantify
fastq_dir = item.metadata["fastq_dir"]
abundance = quantify(fastq_dir)
item.metadata["abundance"] = abundance
manager.mark_done(item)
```

**Tip**: Keep metadata serializable (JSON-compatible) for potential checkpoint persistence.

## Stage Lifecycle

Stages (from `Stage` enum):

| Stage | Meaning | TUI Color |
|-------|---------|-----------|
| `PENDING` | Queued, not yet started | Blue |
| `RUNNING` | Currently processing | Cyan |
| `DONE` | Completed successfully | Green |
| `FAILED` | Terminated with error | Red |

**Backward compatibility**: `SampleStage` (RNA-seq specific) extends this with `DOWNLOADING`, `DOWNLOADED`, `EXTRACTING`, `EXTRACTED`, `QUANTIFYING`, `TRUNCATED`, `SKIPPED`.

**Transition rules**:
- Initial: `PENDING`
- Handler must call `mark_running()` → `RUNNING` (optional but recommended)
- On success: `mark_done()` → `DONE`
- On failure: `mark_failed()` → `FAILED`
- `SKIPPED`/`TRUNCATED` are domain-specific variants of non-success

## Phase Handler Patterns

### Pattern 1: Sequential (Blocking)

Simple loop; each item processed fully before next:

```python
def sequential_handler(manager: BasePipelineManager, items: list[PipelineItem]) -> None:
    for item in items:
        manager.mark_running(item, status="Working...")
        result = process(item.item_id)  # blocking call
        if result.ok:
            manager.mark_done(item, status="Done")
        else:
            manager.mark_failed(item, result.error)
```

**Use for**: Fast operations (< 1s per item), low parallelism needed.

### Pattern 2: Parallel Fan-Out (Non-blocking)

Submit all items to executor, return immediately; completion marks items:

```python
def parallel_handler(manager: BasePipelineManager, items: list[PipelineItem]) -> None:
    futures = {}
    for item in items:
        future = manager.executor.submit(process_one, item)
        futures[future] = item
        manager.mark_running(item, status="Queued")

    # Poll progress (optional)
    while not all(f.done() for f in futures):
        update_progress(manager, futures)
        time.sleep(0.5)

    # Collect results
    for future in as_completed(futures):
        item = futures[future]
        try:
            result = future.result()
            manager.mark_done(item, status="Completed")
        except Exception as exc:
            manager.mark_failed(item, str(exc))
```

**Use for**: I/O-bound (downloads, API calls) or CPU-bound parallelizable work.

**Key**: Handler returns **only after all futures complete** (blocking until all done). If you need truly asynchronous phases, redesign pipeline (BasePipelineManager assumes sequential phase progression).

### Pattern 3: Progress-Polling (Long-Running)

For external processes (subprocess.Popen) or large downloads:

```python
def download_handler(manager: BasePipelineManager, items: list[PipelineItem]) -> None:
    procs = {}
    for item in items:
        cmd = ["wget", item.metadata["url"], "-O", dest]
        proc = subprocess.Popen(cmd, stdout=PIPE, stderr=PIPE)
        procs[proc] = item
        manager.mark_running(item, status="Downloading")

    # Poll all processes
    while procs:
        for proc, item in list(procs.items()):
            retcode = proc.poll()
            if retcode is None:
                # Still running — report bytes downloaded via file size
                current = Path(dest).stat().st_size
                manager.ui.update(item.item_id, current=current / 1e6)
            else:
                # Finished
                del procs[proc]
                if retcode == 0:
                    manager.mark_done(item, status="Downloaded")
                else:
                    manager.mark_failed(item, f"exit {retcode}")
        time.sleep(0.5)
```

**See**: `_rna_download_phase()` in `workflow_manager.py` for production implementation using `robust_download_url()` and file-size polling.

### Pattern 4: Conditional Branching (Dynamic Next Phase)

Rare case where phase order is not fixed:

```python
def quality_check_handler(manager: BasePipelineManager, items: list[PipelineItem]) -> None:
    for item in items:
        qc_passed = run_qc(item.metadata["output"])
        if qc_passed:
            item.metadata["next_phase"] = "analysis"
        else:
            item.metadata["next_phase"] = "reprocess"

def rerun_phase_filter(item: PipelineItem) -> bool:
    # Custom filter used by later phases
    return item.metadata.get("next_phase") == "reprocess"
```

`BasePipelineManager` does not natively support dynamic phase reordering. For dynamic workflows, either:
- Use config-driven workflow (which can define multiple workflows)
- Write a custom manager that adjusts `self.phases` between phases

## Threading and Concurrency

### Executor Model

`BasePipelineManager` owns a single `ThreadPoolExecutor` shared by all phase handlers. All handlers run in the **main thread** (the thread that called `run()`).

**Concurrency achieved by**: Handlers submit tasks to `manager.executor`; tasks run in worker threads; handlers block on `as_completed()` or return immediately (if asynchronous completion callbacks used).

### Thread Safety Rules

1. **PipelineItem access**: An item is only accessed by one handler at a time (phases sequential). However, within a parallel handler:
   - Main thread owns item objects (modifying `stage`, `metadata`, `error`)
   - Worker threads may hold item in closure but should not mutate it

2. **UI updates**: `manager.ui.update()` is thread-safe (implemented with locks). Safe to call from worker threads.

3. **Shared resources**: If accessing shared state (e.g., global counters), protect with `threading.Lock`.

### Avoiding Deadlocks

- Never call `manager.run()` recursively (nested pipeline execution) on same manager instance
- If a handler needs to spawn sub-pipelines, create a **new** `BasePipelineManager` instance:

```python
def parent_handler(manager, items):
    for item in items:
        child = ChildWorkflow(item.metadata["config"])
        child_results = child.run()  # independent executor
        item.metadata["child"] = child_results
        manager.mark_done(item)
```

### Resource-Aware Worker Count

Don't hardcode `max_threads`. Use heuristics:

```python
from metainformant.core.execution.parallel import resource_aware_workers

workers = resource_aware_workers(
    task_type="io",              # "io" or "cpu"
    max_cap=32,                  # Hard upper limit
    memory_per_worker_mb=256,    # Estimated per-thread memory
)
manager = BasePipelineManager(phases, max_threads=workers)
```

**Defaults**:
- `task_type="io"` → `min(cores * 4, 32)`
- `task_type="cpu"` → `max(1, cores - 1)`

Memory constraint (if `psutil` available):
```python
available_mb = psutil.virtual_memory().available / 1e6
workers = min(workers, int(available_mb * 0.7 / memory_per_worker_mb))
```

## TUI Integration

The pipeline manager automatically provides a rich Terminal UI.

### Starting/Stopping (Automatic)

```python
manager.run()  # starts TUI in main thread, blocks until complete
# TUI started before first phase, stopped after last phase (finally block)
```

**Manual control** (rare):

```python
manager.ui.start()
for phase in manager.phases:
    phase.handler(manager, eligible)
manager.ui.stop()
```

### UI Elements

```

 Sample_001 [] 65% 120 MB/s
 Sample_002 [] 100% 45 MB/s
 Sample_003 [] 30% 12 MB/s

 Phase: Download | Active: 3 | Done: 0 | Total: 3 | 65%

```

### Updating Progress

```python
manager.ui.update(
    item_id="Sample_001",
    current=125.5,      # Current bytes/MB/units (depends on context)
    total=1000.0,       # Total expected (for percentage bar)
    status="Downloading",  # Short status text
    color=CYAN,           # Bar color
    speed="1.2 MB/s"      # Optional speed indicator
)
```

**Convenience methods** (via `mark_*`):

```python
manager.mark_running(item, status="Processing...", color=CYAN)
manager.mark_done(item, status="Done", color=GREEN)
manager.mark_failed(item, "Network error", status="Failed", color=RED)
```

These call `ui.update()` internally with sensible defaults (stage → color mapping, status strings).

### Footer Updates

```python
manager.ui.set_footer(f"Phase: {phase.name} | Processed {n}/{total}")
```

Typically called at start of each phase handler.

## Configuration-Driven Workflows

For declarative pipelines, use `metainformant.core.execution.workflow`.

### Validate a Config

```python
from metainformant.core.execution.workflow import validate_config_file

is_valid, errors = validate_config_file("pipeline.yaml")
if not is_valid:
    for error in errors:
        print(f"ERROR: {error}")
```

**Schema**: Supports both `steps:` dict-style (generic) and `downloads:` + `processing:` (two-phase).

### Generate Template

```python
from metainformant.core.execution.workflow import create_sample_config

create_sample_config("my_pipeline.yaml", sample_type="basic")
# Options: "basic", "scientific", "advanced"
```

### Execute Config

```python
from metainformant.core.execution.workflow import download_and_process_data

results = download_and_process_data(
    config_data="pipeline.yaml",  # or dict directly
    output_dir="output/experiment1",
    verbose=True
)

print(f"Downloads: {results['downloads']}")
print(f"Processing: {results['processing']}")
print(f"Errors: {results['errors']}")
```

### Config Schema

```yaml
# pipeline.yaml
name: "My Bioinformatics Pipeline"
description: "Download and analyze expression data"

downloads:
  reference_genome:
    url: "https://example.com/genome.fa.gz"
    filename: "genome.fa.gz"
    checksum: "sha256:abcd1234..."  # optional integrity check

  expression_data:
    url: "https://example.com/expr.tsv"
    filename: "expression.tsv"

processing:
  align:
    type: "bwa_mem"
    reference: "output/download/genome.fa.gz"
    threads: 8

  quantify:
    type: "feature_counts"
    gtf: "annotations.gtf"
```

**Supported sections**:
- `downloads: dict[str, dict]` — Each has `url`, `filename` (optional), `checksum` (optional)
- `processing: dict[str, dict]` — Each has `type` (string) and parameters

### Config Inheritance (Advanced)

Combine base config with overrides programmatically:

```python
from metainformant.core.utils.config import merge_configs

base = load_mapping_from_file("base.yaml")
override = load_mapping_from_file("experiment.yaml")
config = merge_configs(base, override)
```

See `src/metainformant/core/utils/config.py` for merge strategies (deep merge, type coercion).

## Advanced Techniques

### Custom Stage Transitions (Nonstandard Stages)

```python
from metainformant.core.engine.workflow_manager import Stage

# Extend with domain-specific stages
class MyStage(Enum):
    PENDING = "Pending"
    VALIDATING = "Validating"
    VALIDATED = "Validated"
    PROCESSING = "Processing"
    DONE = "Done"
    FAILED = "Failed"

# Use custom stages in handlers
def validate_handler(manager, items):
    for item in items:
        item.stage = MyStage.VALIDATING
        manager.ui.update(item.item_id, stage="Validate", status="Checking...")
        # ... validation ...
        item.stage = MyStage.VALIDATED
        manager.mark_done(item)  # still calls DONE — you could create custom mark method
```

**Better**: Subclass `PipelineItem` and override `mark_done()` to map custom stage→final status.

### Checkpoint / Restart

Persist pipeline state to disk between runs:

```python
import json
from pathlib import Path

def save_checkpoint(manager: BasePipelineManager, path: Path):
    state = {
        item_id: {"stage": item.stage.value, "metadata": item.metadata}
        for item_id, item in manager.items.items()
    }
    io.dump_json(state, path)

def load_checkpoint(path: Path) -> dict:
    return io.load_json(path)

# Resume workflow
checkpoint = load_checkpoint("checkpoint.json")
for item_id, state in checkpoint.items():
    item = PipelineItem(item_id, metadata=state["metadata"])
    item.stage = Stage(state["stage"])
    manager.items[item_id] = item
```

**Caveat**: `PipelinePhase.filter_fn` must recognize `Stage.DONE` as already finished (default filter only matches `PENDING`).

### Dynamic Phase Addition (Runtime)

Modify `manager.phases` before calling `run()`:

```python
# Conditional phase inclusion
if config.get("enable_secondary_analysis"):
    manager.phases.append(PipelinePhase("Secondary", secondary_handler))
```

Cannot add phases mid-run without careful state management (not recommended for typical use).

### Handler Composition (Nested Calls)

Break complex handlers into sub-functions:

```python
def phase_handler(manager, items):
    for item in items:
        manager.mark_running(item)
        # Composition: smaller unit-testable functions
        preprocess(item)
        transformed = transform(item)
        postprocess(transformed)
        manager.mark_done(item)
```

Keep functions pure where possible; side effects limited to `manager` and `item.metadata`.

## Common Pitfalls

### Pitfall 1: Handler Never Calls `mark_done()` / `mark_failed()`

**Symptom**: Pipeline hangs indefinitely after phase output stops updating.

**Cause**: Handler return without marking item completion; `run()` waits for all eligible items to finish.

**Fix**: Ensure every code path in handler calls `mark_done()` or `mark_failed()` for each item:

```python
def handler(manager, items):
    for item in items:
        try:
            do_work(item)
            manager.mark_done(item)
        except Exception as exc:
            manager.mark_failed(item, str(exc))  # catches and marks
```

### Pitfall 2: Exceptions Thrown in Worker Threads Not Propagated

**Symptom**: Pipeline completes but some items failed silently; `results[item_id] == False`.

**Cause**: Worker thread exceptions captured by `Future`; main thread never inspects them.

**Fix**: In handler, always retrieve `future.result()` inside try/except:

```python
futures = {manager.executor.submit(fn, item): item for item in items}
for future in as_completed(futures):
    item = futures[future]
    try:
        future.result()  # re-raises exception
        manager.mark_done(item)
    except Exception as exc:
        manager.mark_failed(item, str(exc))
```

### Pitfall 3: Modifying `item.stage` Directly Without `mark_*`

**Symptom**: TUI doesn't update, `results` dict incorrect.

**Cause**: Direct assignment like `item.stage = Stage.DONE` bypasses TUI updates and bookkeeping.

**Fix**: Always use manager methods:

```python
# CORRECT
manager.mark_done(item)

# WRONG
item.stage = Stage.DONE  # TUI won't reflect this
```

### Pitfall 4: Too Many Threads Causing OOM

**Symptom**: System slows to a halt, process killed by OOM killer.

**Cause**: `max_threads` too high for memory-intensive tasks; each thread loads large dataset.

**Fix**:
- Use `resource_aware_workers(memory_per_worker_mb=...)`
- Or divide items into batches processed sequentially

```python
# Batch processing: 100 items, 5 workers, 20 batches
for batch in chunked(items, 20):
    manager = BasePipelineManager(phases, max_threads=5)
    for item in batch:
        manager.add_item(item.id)
    manager.run()
```

### Pitfall 5: Blocking I/O in Main Thread

**Symptom**: TUI freezes during phase; progress bar stops updating.

**Cause**: Handler performs network I/O or subprocess.wait() on main thread.

**Fix**: Offload blocking work to `manager.executor`:

```python
def handler(manager, items):
    futures = [manager.executor.submit(blocking_io, item) for item in items]
    for future in as_completed(futures):
        # This loop yields to allow TUI to refresh
        result = future.result()
        manager.mark_done(...)
```

## Testing Orchestration Code

### Test Philosophy (Real Implementation)

All tests use real implementations:
- Real file I/O to `tmp_path` (fixture)
- Real subprocesses (no `test-double patching APIs` on `subprocess.Popen`)
- Real thread pools

See: [../../tests/REAL_IMPLEMENTATION_TESTING_POLICY.md](../../tests/REAL_IMPLEMENTATION_TESTING_POLICY.md)

### Example Test

```python
def test_simple_pipeline(tmp_path):
    # Arrange
    phases = [
        PipelinePhase("Write", write_handler),
        PipelinePhase("Read", read_handler),
    ]
    manager = BasePipelineManager(phases, max_threads=2)
    manager.add_item("test", metadata={"output": tmp_path / "out.txt"})

    # Act
    results = manager.run()

    # Assert
    assert results["test"] is True
    assert (tmp_path / "out.txt").exists()
```

**Test-specific TUI disable**: `TerminalInterface` can run in headless mode; set env var `METAINFORMANT_NO_TUI=1` or pass `use_tui=False` to manager (if implemented).

### Integration Test with Real External Commands

```python
def test_download_real_sra(tmp_path):
    phases = [PipelinePhase("Download", _rna_download_phase)]  # uses robust_download_url
    manager = BasePipelineManager(phases, max_threads=1)
    manager.add_item("ERR123456", metadata={"url": "https://..."})
    results = manager.run()
    # Actually downloads! (tests should use small public datasets)
```

## Comparison: `BasePipelineManager` vs Config-Driven

| Feature | BasePipelineManager | Config-driven (workflow.py) |
|---------|--------------------|-----------------------------|
| API style | OOP, Python code | Declarative YAML/JSON |
| Flexibility | Full programmability | Limited to download/process schema |
| User friendliness | Developer-facing | End-user-facing |
| Conditional logic | Yes (in handler) | No (static) |
| Phase order | Fixed at construction | Fixed by schema (downloads then processing) |
| Parallelism | Explicit via `executor.submit()` | Implicit in `download_and_process_data` (thread pool internally) |
| Phase count | Arbitrary | Implicitly 2 (download, processing) |
| Extensibility | Custom phases, filters, TUI | Need code modifications |

**Guideline**:
- For new domain modules with >2 steps → Use `BasePipelineManager`
- For simple 2-step (fetch then compute) → Config-driven is quicker

## Future Directions

- **Async/Await**: Potential migration to `asyncio` for better I/O concurrency
- **Distributed Execution**: Workers across multiple nodes (Redis/Celery backend)
- **DAG Workflows**: Phase dependencies as DAG instead of linear sequence
- **Checkpoint Persistence**: Automatic state serialization to resume interrupted runs
- **Dynamic Scaling**: Elastic worker pool adjusting to load

## See Also

- [Architecture](ARCHITECTURE.md) — System context and coordination patterns
- [Multi-Agent Workflows](MULTI_AGENT_WORKFLOWS.md) — Real-world pipeline compositions
- [Communication Protocols](COMMUNICATION_PROTOCOLS.md) — Data flow between agents
- [Safety](SAFETY.md) — Error handling and validation
- Source code: [`src/metainformant/core/engine/workflow_manager.py`](../../src/metainformant/core/engine/workflow_manager.py)
- Source code: [`src/metainformant/core/execution/workflow.py`](../../src/metainformant/core/execution/workflow.py)
- Source code: [`src/metainformant/core/execution/parallel.py`](../../src/metainformant/core/execution/parallel.py)
