# Troubleshooting: Agent Coordination Issues

**Troubleshooting guide for multi-agent workflow failures and coordination problems**

## Quick Diagnosis by Symptom

| Symptom | Likely Cause | First Check | Fix |
|---------|--------------|-------------|-----|
| Pipeline hangs indefinitely | Handler never calls `mark_done()` | Is any item stuck in `RUNNING` in `checkpoint.json`? | Add timeout or ensure all code paths call `mark_done()` |
| Some items fail, others succeed | Per-item error isolated | `item.error` in checkpoint or logs | Inspect error; fix input or skip bad items |
| Memory exhaustion | Too many parallel workers | `ps aux | grep python` memory usage | Reduce `max_threads`; use `resource_aware_workers()` |
| Disk fills up | Unchecked intermediate files | `df -h output/` | Enable compression; prune old runs |
| TUI freezes | Blocking I/O on main thread | `top` shows Python at 100% CPU but no progress | Offload blocking ops to `manager.executor.submit()` |
| Items fail with "FileNotFound" | Downstream phase runs before upstream completes | Check if file exists in expected `output/` path | Verify phase ordering; check `item.metadata["path"]` |
| All items fail immediately | Systemic error (missing dependency, bad config) | First error in logs | Check config, external tool installation |
| Slow progress despite many workers | Contention (disk I/O, network bandwidth) | `iotop` shows saturation | Reduce parallelism or stagger phases |
| Results inconsistent between runs | Non-deterministic code (random seeds) | Compare output hashes | Set `random.seed()` / `np.random.seed()` |
| Restart doesn't resume | Checkpoint file corrupted or partial | Inspect `checkpoint.json` JSON syntax | Remove checkpoint and restart from scratch |
| Missing metadata keys in phase B | Phase A didn't set key or early exit | Log `item.metadata` at end of Phase A | Ensure Phase A sets all expected keys even on partial success |

## Detailed Troubleshooting

### 1. Pipeline Stalled / No Progress

**Symptom**: TUI active but progress bars frozen for > 30 seconds.

**Diagnostic**:

```bash
# 1. Check which items are RUNNING but not DONE
python -c "
import json
state = json.load(open('output/checkpoint.json'))
for item_id, data in state['items'].items():
    if data['stage'] == 'RUNNING':
        print(f'Stuck: {item_id}')
"

# 2. Tail logs for warnings around stall time
tail -f output/core/logs/pipeline.log | grep -i "stuck\|timeout\|blocking"
```

**Common causes**:
- Handler submitted futures but never waited for `as_completed()` → items stuck RUNNING forever
- Handler blocked on network I/O on main thread (no TUI updates)
- Infinite loop or long computation without yielding

**Fix**:
- Ensure parallel handlers call `future.result()` inside try/except and handle `TimeoutError`
- Add `time.sleep(0.5)` in polling loops to yield
- Move blocking I/O to executor threads

**Code template**:

```python
def handler(manager, items):
    futures = [manager.executor.submit(task, item) for item in items]
    for future in as_completed(futures, timeout=300):
        try:
            future.result()  # will raise if task failed
            manager.mark_done(...)
        except TimeoutError:
            manager.mark_failed(..., "TIMEOUT")
```

---

### 2. Memory Leak / OOM Killer

**Symptom**: System slows, swap activity high, process killed.

**Diagnostic**:

```bash
# Monitor memory usage over time
watch -n 1 'ps -o pid,rss,cmd -p $(pgrep -f metainformant)'

# Check kernel messages for OOM kill
dmesg | grep -i "killed process"
```

**Common causes**:
- Loading all items into memory at once (`io.load_json()` huge file)
- Metadata storing large objects (DataFrames)
- Thread pool with thousands of workers (threads = items)

**Fix**:
- Stream large files (line-by-line for JSONL)
- Delete objects explicitly: `del large_obj; import gc; gc.collect()`
- Use `resource_aware_workers()` to cap threads based on available memory
- Store large results on disk, pass paths in metadata not objects

---

### 3. Cascading Failures (All Items Fail)

**Symptom**: Every item ends as `FAILED` in same phase.

**Diagnostic**:

```bash
# Look at log for phase start and first error
grep "Phase: X" output/core/logs/pipeline.log | head -1
grep -A 5 "ERROR.*Phase X" output/core/logs/pipeline.log | head -20
```

**Common causes**:
- External dependency not installed or not in PATH
- Config file missing required key
- Permission denied on output directory
- Network unreachable for all downloads

**Fix**:
- Validate external tool: `which kallisto` returns path
- Validate config: `python -m metainformant.core.execution.workflow --validate config.yaml`
- Check write permissions: `touch output/test.txt` (no error)
- Ping external host: `ping ncbi.nlm.nih.gov`

---

### 4. Intermittent Network Errors (Retry Exhausted)

**Symptom**: Items fail with connection errors after N retries.

**Diagnostic**:

```bash
# Sample error messages
grep "TIMEOUT\|ConnectionError\|HTTP" output/core/logs/*.log
```

**Common causes**:
- Flaky network (WiFi dropouts)
- Remote server rate-limiting
- Firewall blocking

**Fix**:
- Increase `retry_attempts` in config
- Add exponential backoff (already in place for downloads)
- Use ENA mirrors if NCBI slow
- Run on stable network (Ethernet, not WiFi)

---

### 5. Output Missing or Incomplete

**Symptom**: Pipeline completes (DONE) but expected output files not present.

**Diagnostic**:

```bash
# Search output for item ID
find output/ -name "*SRR123*" -type f

# Check metadata for that item
grep -A 10 '"SRR123"' output/checkpoint.json
```

**Common causes**:
- Handler marked DONE without actually writing file
- File written to wrong directory (relative path vs absolute)
- Atomic write failed silently (disk full)

**Fix**:
- Add verification step before `mark_done()`:

  ```python
  if not output_path.exists():
      raise FileNotFoundError(f"Expected output missing: {output_path}")
  manager.mark_done(item)
  ```

- Use absolute paths: `output_path = Path(manager.config["output_dir"]) / item_id / "file.json"`
- Check disk space before phase: `ensure_disk_space(needed_gb=10)`

---

### 6. Metadata Key Not Found (KeyError in downstream phase)

**Symptom**: Phase B handler crashes with `KeyError: 'some_key'`.

**Diagnostic**:

```python
# Add debug print at start of Phase B
def phase_b(manager, items):
    for item in items:
        print(f"DEBUG: metadata keys = {list(item.metadata.keys())}")
        # ...
```

**Common causes**:
- Phase A didn't set the key (typo in key name)
- Phase A failed but Phase B ran (filter didn't exclude FAILED items)
- Handler mutated metadata (deleted key)

**Fix**:
- Use consistent key names; consider namespace: `"rna:fastq_dir"` vs `"fastq"` across modules
- Ensure default filter (`lambda i: i.stage == Stage.PENDING`) excludes non-PENDING items
- Add defensive code in Phase B:

  ```python
  if "required_key" not in item.metadata:
      manager.mark_failed(item, "Missing required_key from previous phase")
      continue
  ```

---

### 7. TUI Characters Garbled

**Symptom**: Terminal shows weird characters, progress bars misaligned.

**Cause**: stdout is not a TTY (output piped to file or CI environment).

**Fix**:

```bash
# Disable TUI
export METAINFORMANT_NO_TUI=1

# Or programmatically
manager = BasePipelineManager(phases, use_tui=False)
```

---

### 8. Progress Not Updating (TUI frozen but handler running)

**Symptom**: TUI shows last update; main thread appears stuck but actually working.

**Diagnostic**:

```bash
# In another terminal, find Python process and check threads
ps -T -p <pid> | grep -i thread
```

**Cause**: TUI refresh loop runs on main thread; if main thread blocked in long computation, UI cannot redraw.

**Fix**: Offload blocking work to executor:

```python
# WRONG: blocking call on main thread
result = subprocess.run(cmd, check=True)  # blocks UI

# RIGHT: submit to executor
future = manager.executor.submit(subprocess.run, cmd, check=True)
# Main thread can continue TUI loop; callback marks done later
```

---

### 9. Checkpoint Restore Doesn't Work

**Symptom**: After restore, pipeline still re-runs completed items.

**Diagnostic**:

```python
# Inspect checkpoint structure
import json
ckpt = json.load(open("checkpoint.json"))
print(ckpt["items"]["SRR123"]["stage"])  # should be "DONE"
```

**Cause**: Phase filter doesn't recognize restored `Stage.DONE` as already completed.

**Fix**: Ensure phase `filter_fn` excludes DONE items:

```python
# Default: only PENDING items are eligible
PipelinePhase("MyPhase", handler)  # filter_fn = lambda i: i.stage == Stage.PENDING

# Custom: include already-DONE to skip (but better to remove phase entirely)
PipelinePhase("MyPhase", handler, filter_fn=lambda i: i.stage != Stage.FAILED)
```

**Better**: Filter out DONE items from `manager.phases` before calling `run()`.

---

### 10. Slow Performance Despite High Worker Count

**Symptom**: `max_threads=32` but only ~5% CPU used.

**Diagnostic**:

```bash
# Check I/O wait
iostat -x 1

# Check network
iftop -n

# Check disk queue
iostat -d -x 1
```

**Cause**: Bottleneck outside CPU (disk I/O, network bandwidth, external process).

**Fix**:
- **Disk I/O bound**: Reduce parallelism; process items in batches
- **Network bound**: Already maxed out by few threads; more workers won't help
- **External tool single-threaded**: Don't parallelize beyond tool's capacity (e.g., `kallisto quant` is single-threaded per sample)

Increase resources: faster SSD, higher bandwidth network, tool with threading flag.

---

## Debug Commands Cheat Sheet

```bash
# 1. View aggregated logs
tail -f output/core/logs/*.log | grep -E "ERROR|CRITICAL|WARNING"

# 2. Inspect live checkpoint
jq '.items | to_entries | map(select(.value.stage=="FAILED"))' output/checkpoint.json

# 3. Count items by stage
jq '[.items[].stage] | group_by(.) | map({stage: .[0], count: length})' output/checkpoint.json

# 4. Check disk usage
du -sh output/* | sort -hr  # largest output dirs

# 5. Monitor memory
watch -n 1 'ps -o pid,rss,cmd -p $(pgrep -f metainformant)'

# 6. Validate a config file
python -c "
from metainformant.core.execution.workflow import validate_config_file
ok, errs = validate_config_file('config.yaml')
print('OK' if ok else '\\n'.join(errs))
"

# 7. Dry-run (if supported by phase)
python -m metainformant.rna.workflow --dry-run --config config.yaml
```

---

## When to Escalate

| Condition | Escalate to | Provide |
|-----------|-------------|---------|
| External tool consistently crashes | Domain module maintainer | Tool version, command, full stderr |
| Workflow deadlocks across all modules | Core maintainer (`core.engine`) | Checkpoint JSON, logs, `htop` screenshot |
| TUI crashes with traceback | UI maintainer (`core.ui`) | Terminal type (`echo $TERM`), OS, Python version |
| Checkpoint restore produces wrong results | Workflow maintainer | Pre- and post-checkpoint item states |
| Performance outlier (10× slower than baseline) | Performance team | Config, resource stats (`vmstat 1`), thread count |

---

## Common Error Messages and Resolutions

| Error Message | Meaning | Likely Fix |
|---------------|---------|-----------|
| `OSError: [Errno 28] No space left on device` | Disk full | Clean `output/`; mount larger volume |
| `NetworkError: Connection refused` | Remote host unreachable | Check firewall; retry later |
| `ValueError: Checksum mismatch` | Corrupted or tampered file | Delete and re-download; check mirror |
| `ImportError: cannot import name 'X'` | Incompatible version | `uv sync --refresh` to update |
| `PermissionError: [Errno 13]` | File not writable | `chmod u+w output/` or choose different path |
| `subprocess.CalledProcessError: command '...' returned non-zero` | External tool failed | Inspect command stderr in logs for tool-specific message |
| `TimeoutError` | Operation took > configured timeout | Increase timeout in config; check system load |
| `json.JSONDecodeError` | Invalid JSON file | Inspect file; may be partial write — consider atomic write issue |
| `FileNotFoundError` for expected output | Upstream phase didn't produce file | Check upstream logs; verify path passed in metadata is correct |

---

## Checklist for Debug Sessions

- [ ] Enable DEBUG log level: `export METAINFORMANT_LOG_LEVEL=DEBUG`
- [ ] Inspect `output/checkpoint.json` for partial state
- [ ] Tail logs: `tail -f output/core/logs/pipeline.log`
- [ ] Validate config schema
- [ ] Verify external tools in PATH: `which kallisto plink bwa`
- [ ] Check disk space: `df -h output/`
- [ ] Check memory: `free -h`
- [ ] Run workflow with single item to reproduce in isolation
- [ ] Search logs for first ERROR line (root cause often there)
- [ ] Re-run with `METAINFORMANT_NO_TUI=1` to eliminate TUI as suspect

---

## Next Steps

After resolving the issue:
1. Document root cause in code comment or issue tracker
2. Add validation or guard if preventable (e.g., disk space check)
3. Consider adding unit test that reproduces scenario
4. Update [Safety](SAFETY.md) or [Best Practices](BEST_PRACTICES.md) if pattern discovered

---

*Happy debugging! Remember: Every failure is an opportunity to make the system more robust.*
