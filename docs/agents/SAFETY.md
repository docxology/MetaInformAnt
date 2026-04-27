# Safety: Error Handling, Validation, and Recovery

**Status**: Implementation guide  
**Audience**: All developers building agents or workflows  
**Scope**: Protecting data integrity, preventing cascading failures, enabling recovery

## Table of Contents

- [Philosophy](#philosophy)
- [Validation Layers](#validation-layers)
- [Error Isolation](#error-isolation)
- [Atomic Operations](#atomic-operations)
- [Retry Policies](#retry-policies)
- [Recovery and Rollback](#recovery-and-rollback)
- [Resource Limits](#resource-limits)
- [Monitoring and Alerts](#monitoring-and-alerts)
- [Checkpoint/Restart](#checkpointrestart)
- [Best Practices Checklist](#best-practices-checklist)

## Philosophy

METAINFORMANT adopts a **fail-safe, fail-separable** philosophy:

1. **Fail safe**: Errors should not leave partial/invalid outputs
2. **Fail separable**: One agent's failure should not cascade to others
3. **Diagnosable**: Errors must be logged with enough context to reproduce
4. **Recoverable**: Where possible, workflows can resume from last good state

### Principles in Code

```python
def safe_handler(manager: BasePipelineManager, items: list[PipelineItem]) -> None:
    for item in items:
        manager.mark_running(item)
        try:
            # Validate inputs
            validate_inputs(item)

            # Do work atomically
            result = atomic_operation(item)

            # Verify output
            if not verify_output(result):
                raise ValueError("Output validation failed")

            # Write to temp, then atomic rename
            temp_path = output_path.with_suffix(".tmp")
            io.dump_json(result, temp_path)
            temp_path.rename(output_path)  # atomic

            manager.mark_done(item)

        except (ValueError, IOError) as exc:
            logger.error(
                "Item failed",
                extra={"item_id": item.item_id, "error": str(exc)},
                exc_info=True,
            )
            manager.mark_failed(item, f"SAFE_ERROR: {exc}")
```

---

## Validation Layers

Validation occurs at multiple boundaries to catch errors early.

### Layer 1: Input Validation (Gatekeeper)

Before any computation, validate inputs:

```python
def validate_inputs(item: PipelineItem) -> None:
    """Check that required metadata keys exist and are well-formed."""
    required = ["file_path", "sample_id"]
    for key in required:
        if key not in item.metadata:
            raise ValueError(f"Missing required key: {key}")

    path = Path(item.metadata["file_path"])
    if not path.exists():
        raise FileNotFoundError(f"Input file missing: {path}")

    if path.stat().st_size == 0:
        raise ValueError("Input file is empty")

    # Optional: checksum verification
    expected = item.metadata.get("expected_checksum")
    if expected:
        actual = io.checksums.sha256sum(path)
        if actual != expected:
            raise ValueError(f"Checksum mismatch: expected {expected}, got {actual}")
```

**Use**: First step in every phase handler.

### Layer 2: Config Validation

```python
from metainformant.core.execution.workflow import validate_config_file

is_valid, errors = validate_config_file("pipeline.yaml")
if not is_valid:
    for error in errors:
        logger.critical(f"Config invalid: {error}")
    sys.exit(1)
```

**Schema enforcement**: JSON Schema validation (if `schema_path` provided).

### Layer 3: Runtime Assertions

```python
def compute(item):
    result = heavy_computation()
    assert result is not None, "Computation returned None"
    assert result > 0, "Result must be positive"
    return result
```

Assertions fail fast during development; remove or convert to exceptions in production if needed.

### Layer 4: Output Verification

After computation, verify outputs before marking DONE:

```python
def verify_output(output_path: Path) -> bool:
    """Check file exists, non-empty, valid structure."""
    if not output_path.exists():
        return False
    if output_path.stat().st_size == 0:
        return False
    # Optional: parse to validate format
    try:
        _ = io.load_json(output_path)
        return True
    except Exception:
        return False
```

---

## Error Isolation

### Per-Item Error Boundary

Each `PipelineItem` is an independent error domain.

#### Implementation

```python
def handler(manager, items):
    for item in items:
        try:
            process(item)
            manager.mark_done(item)
        except Exception as exc:
            # Item fails in isolation; other items continue
            manager.mark_failed(item, str(exc))
```

**Result**: Results dict:

```python
{
    "item1": True,   # success
    "item2": False,  # failed; error in item.error
    "item3": True,
}
```

### Per-Phase Error Aggregation

If a phase handler spawns parallel workers:

```python
def parallel_handler(manager, items):
    futures = {manager.executor.submit(fn, item): item for item in items}
    errors = []
    for future in as_completed(futures):
        item = futures[future]
        try:
            future.result()  # re-raise exception from worker
            manager.mark_done(item)
        except Exception as exc:
            error_msg = f"{item.item_id}: {exc}"
            errors.append(error_msg)
            manager.mark_failed(item, str(exc))

    # Log summary at phase end
    if errors:
        logger.error(f"Phase completed with {len(errors)} failures")
```

### Error Propagation vs. Suppression

| Scenario | Action |
|----------|--------|
| Single item failure | Isolate; mark FAILED; continue |
| All items fail | Consider aborting pipeline early (after phase) |
| External dependency down | Fail fast: mark all items FAILED, don't retry indefinitely |
| Transient network error | Retry with backoff (see Retry Policies) |

---

## Atomic Operations

### Atomic File Writes

Use `io.dump_json()` with `atomic=True` (default):

```python
io.dump_json(data, "output/result.json")  # writes to .tmp then renames
```

Implementation:

```python
# Pseudocode of atomic=True behavior
tmp_path = output_path.with_suffix(".tmp")
with open(tmp_path, "w") as f:
    json.dump(data, f)
tmp_path.rename(output_path)  # atomic on POSIX
```

### Checkpoint Atomicity

Write checkpoint in one shot:

```python
checkpoint = {item_id: item.stage.value for item_id, item in manager.items.items()}
io.dump_json(checkpoint, "checkpoint.json")  # atomic
```

### Database Transactions (if using PostgreSQL)

When `metainformant.core.data.db` is used, wrap mutations in transactions:

```python
from metainformant.core.data import db

with db.transaction():
    db.insert_results(results)
    db.log_phase_completion(phase_name, item_count)
# Commit on success; rollback on exception
```

---

## Retry Policies

### Exponential Backoff for Transient Errors

```python
import time
from tenacity import retry, stop_after_attempt, wait_exponential  # optional dependency

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def download_with_retry(url: str, dest: Path) -> None:
    io.download_file(url, dest)
```

**If `tenacity` unavailable**:

```python
def download_with_retry(url, dest, max_attempts=3):
    for attempt in range(1, max_attempts + 1):
        try:
            io.download_file(url, dest)
            return
        except (NetworkError, Timeout) as exc:
            if attempt == max_attempts:
                raise
            backoff = min(2 ** attempt, 60)  # 2, 4, 8... up to 60s
            time.sleep(backoff)
```

### Per-Item Retry vs. Pipeline Retry

| Strategy | When to Use | Implementation |
|----------|-------------|----------------|
| **Per-item retry** | Network I/O, external API calls | Retry inside handler before marking failed |
| **Pipeline-level retry** | Entire phase failed due to systemic issue | Rerun `manager.run()` from checkpoint |
| **No retry** | Irrecoverable errors (syntax, validation) | Immediate fail |

**Recommendation**: Retry transient failures at handler level (max 3 attempts); don't retry permanent errors.

### Circuit Breaker Pattern

For external services (API, database):

```python
class CircuitBreaker:
    def __init__(self, max_failures=5, reset_timeout=60):
        self.failures = 0
        self.last_failure = 0
        self.max_failures = max_failures
        self.reset_timeout = reset_timeout

    def call(self, func, *args, **kwargs):
        if self.failures >= self.max_failures:
            if time.time() - self.last_failure < self.reset_timeout:
                raise CircuitOpenError("Circuit open — service down")
            else:
                self.failures = 0  # half-open: allow one trial

        try:
            result = func(*args, **kwargs)
            self.failures = 0  # success resets
            return result
        except Exception as exc:
            self.failures += 1
            self.last_failure = time.time()
            raise
```

Use for repeated API calls (UniProt, NCBI) to avoid hammering failing endpoint.

---

## Recovery and Rollback

### Strategy 1: No-Op on Duplicate Run (Idempotent)

Design handlers to be idempotent:

```python
def download_handler(manager, items):
    for item in items:
        dest = Path(item.metadata["output_path"])
        if dest.exists() and dest.stat().st_size > 0:
            logger.info(f"{item.item_id} already downloaded — skipping")
            manager.mark_done(item, status="Cached")
            continue

        # Download to temp, then atomic rename
        tmp = dest.with_suffix(".part")
        io.download_file(url, tmp)
        tmp.rename(dest)
        manager.mark_done(item)
```

**Benefit**: Re-running workflow picks up where it left off.

### Strategy 2: Checkpoint Files

Periodic snapshots of `manager.items`:

```python
def checkpoint_handler(manager, items):
    state = {
        item_id: {
            "stage": item.stage.value,
            "metadata": {k: v for k, v in item.metadata.items() if k != "large_object"},
        }
        for item_id, item in manager.items.items()
    }
    io.dump_json(state, "checkpoint.json")
    manager.mark_done(items[0])
```

**Restore**:

```python
state = io.load_json("checkpoint.json")
for item_id, saved in state.items():
    item = manager.items[item_id]
    item.stage = Stage(saved["stage"])
    item.metadata.update(saved.get("metadata", {}))
```

Then resume from phase `i` where stage < DONE.

### Strategy 3: Compensation Actions (Sagas)

For multi-phase distributed transactions, implement undo operations:

```python
class TwoPhaseCommitWorkflow:
    def run(self):
        try:
            self._phase1()  # prepare
            self._phase2()  # commit
        except Exception:
            self._compensate_phase1()  # rollback
            raise
```

**Avoid for single-pipeline**: Prefer idempotency and checkpoint restart over complex sagas.

### Strategy 4: Output Directory Versioning

Tag outputs with run ID (timestamp or UUID):

```python
import uuid
run_id = uuid.uuid4().hex[:8]
output_dir = Path(f"output/runs/{run_id}")
```

Then old outputs preserved; new run creates fresh tree. Crash-safe: old outputs untouched.

---

## Resource Limits

### Memory Guardrails

```python
def memory_guarded_handler(manager, items):
    import psutil
    available = psutil.virtual_memory().available
    if available < 1e9:  # < 1 GB
        logger.warning("Low memory — reducing concurrency")
        manager.max_threads = max(1, manager.max_threads // 2)
```

**Alternatively**: Set process-level limit via `resource` module (Unix):

```python
import resource

# Limit virtual memory to 64 GB
resource.setrlimit(resource.RLIMIT_AS, (64 * 1024**3, 64 * 1024**3))
```

### Disk Space Check

```python
def ensure_disk_space(needed_gb: float, path: Path = Path(".")) -> bool:
    stat = shutil.disk_usage(path)
    free_gb = stat.free / 1e9
    if free_gb < needed_gb * 1.2:  # 20% safety margin
        raise OSError(f"Insufficient disk space: need {needed_gb} GB, have {free_gb:.1f} GB")
    return True
```

Call before large downloads or intensive processing.

### Thread Count Capping

Always compute workers via `resource_aware_workers()`:

```python
workers = resource_aware_workers(task_type="io")
if workers > 32:
    logger.warning(f"Capping workers from {workers} to 32 (safety limit)")
    workers = 32
```

### Timeouts

```python
def handler_with_timeout(manager, items):
    for item in items:
        future = manager.executor.submit(process_one, item)
        try:
            result = future.result(timeout=300)  # 5 minutes max
        except concurrent.futures.TimeoutError:
            manager.mark_failed(item, "TIMEOUT after 300s")
```

**Per-item timeout** prevents one slow item from blocking phase completion.

---

## Monitoring and Alerts

### Structured Logging

All agents should use `get_logger(__name__)` and emit structured logs:

```python
logger.info(
    "Phase completed",
    extra={
        "phase": "download",
        "item_id": item_id,
        "duration_s": duration,
        "bytes": bytes_downloaded,
        "status": "success",
    }
)
```

**Log aggregation**: Capture logs to centralized system (ELK, Loki) for analysis.

### Metrics (Optional)

If `prometheus_client` available:

```python
from prometheus_client import Counter, Histogram

PHASE_DURATION = Histogram("phase_duration_seconds", "Time per phase", ["phase"])
PHASE_FAILURES = Counter("phase_failures_total", "Failures per phase", ["phase"])

def handler(manager, items):
    with PHASE_DURATION.labels(phase="download").time():
        try:
            # work
            pass
        except Exception:
            PHASE_FAILURES.labels(phase="download").inc()
            raise
```

Expose `/metrics` endpoint for monitoring.

### Alert Conditions

Configure alerts on:

- `pipeline_failed_total > 0` (any item failure)
- `phase_duration_seconds{phase="download"} > 3600` (stuck download)
- `disk_free_bytes < 10%` (low disk space)
- `memory_usage_bytes > 90%` (memory pressure)

---

## Checkpoint/Restart

### Checkpoint File Format

```json
{
  "run_id": "20250426_abc123",
  "pipeline": "rna_amalgkit",
  "phases_completed": ["download", "getfastq"],
  "items": {
    "SRR123": {
      "stage": "DONE",
      "metadata": {
        "sra_path": "output/rna/download/SRR123.sra",
        "fastq_dir": "output/rna/fastq/SRR123"
      }
    },
    "SRR124": {
      "stage": "FAILED",
      "error": "Network timeout"
    }
  }
}
```

### Save Checkpoint

```python
def save_checkpoint(manager: BasePipelineManager, path: Path) -> None:
    state = {
        "run_id": manager.config.get("run_id", "unknown"),
        "phases_completed": [p.name for p in manager.phases if ...],  # need tracking
        "items": {
            iid: {
                "stage": item.stage.value,
                "metadata": item.metadata,
            }
            for iid, item in manager.items.items()
        },
    }
    io.dump_json(state, path)
```

### Restart from Checkpoint

```python
def restore_from_checkpoint(
    manager: BasePipelineManager,
    checkpoint_path: Path
) -> None:
    state = io.load_json(checkpoint_path)

    # Rebuild items
    for item_id, item_data in state["items"].items():
        item = PipelineItem(item_id)
        item.stage = Stage(item_data["stage"])
        item.metadata.update(item_data.get("metadata", {}))
        manager.items[item_id] = item

    # Skip phases already completed
    completed_phases = state.get("phases_completed", [])
    manager.phases = [
        phase for phase in manager.phases
        if phase.name not in completed_phases
    ]
```

**Limitation**: `BasePipelineManager` doesn't natively persist phase history; checkpointing requires enhancing manager or storing separately.

### Automatic Checkpointing

Modify `BasePipelineManager.run()` to checkpoint after each phase:

```python
def run_with_checkpoints(self, checkpoint_dir: Path) -> Dict[str, bool]:
    for i, phase in enumerate(self.phases):
        # Run phase
        eligible = [item for item in self.items.values() if phase.filter_fn(item)]
        phase.handler(self, eligible)

        # Checkpoint
        ckpt = checkpoint_dir / f"phase_{i}_{phase.name}.json"
        save_checkpoint(self, ckpt)

    return {iid: item.stage == Stage.DONE for iid, item in self.items.items()}
```

---

## Best Practices Checklist

Before marking a pipeline production-ready, verify:

- [ ] **Input validation** present in all phase handlers
- [ ] **Output verification** before `mark_done()`
- [ ] **Atomic writes** (via `io.dump_json(atomic=True)`)
- [ ] **Resource limits** enforced (`resource_aware_workers` or similar)
- [ ] **Per-item try/except** with `mark_failed()`
- [ ] **Structured logging** at DEBUG level per item
- [ ] **Idempotency**: re-running doesn't duplicate/corrupt outputs
- [ ] **Checkpoint capability** (optional but recommended for long runs)
- [ ] **TUI integration** for progress visibility
- [ ] **Error messages** include actionable context (item ID, file path, error code)

---

## Troubleshooting Safety Failures

| Failure Mode | Diagnosis | Remedy |
|--------------|-----------|--------|
| Partial output file on disk | Handler crashed before atomic rename | Ensure atomic writes (temp→final rename) |
| Disk fills mid-run | `OSError: no space left on device` | Pre-run `ensure_disk_space()`, monitor with alerts |
| Worker thread silent exit | Unhandled exception, future.result() not checked | Always call `future.result()` inside try/except |
| Memory keeps growing | Metadata storing large data, thread pool caching | Store paths not objects; use `executor.shutdown(wait=True)` |
| Stalled pipeline with no logs | Handler blocked on I/O without timeouts | Add timeouts; move I/O to executor threads |
| Cascade failure (all items fail) | Shared dependency down (database, API) | Circuit breaker; batch items into retry groups |

### Debug Session Template

```bash
# 1. Check logs for first failure
grep "ERROR" output/core/logs/*.log | head

# 2. Inspect checkpoint state
cat output/checkpoint.json | jq '.items | map(select(.stage == "FAILED"))'

# 3. Run with debug logging
METAINFORMANT_LOG_LEVEL=DEBUG python -m metainformant.rna.workflow ...

# 4. Validate output files
find output/ -type f -size 0  # empty files?
find output/ -type f -name "*.part"  # incomplete downloads?
```

---

## Next Steps

- Study [Communication Protocols](COMMUNICATION_PROTOCOLS.md) for metadata sharing patterns
- Review [Multi-Agent Workflows](MULTI_AGENT_WORKFLOWS.md) for safety in complex compositions
- Consult [Orchestration](ORCHESTRATION.md) for phase handler API details
- Implement module-specific safety checks in your domain's `rules/{module}.md`
