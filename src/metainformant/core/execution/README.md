# Execution

Parallel execution primitives, config-driven workflow orchestration, and codebase discovery (function indexing, call graphs, symbol search).

## Contents

| File | Purpose |
|------|---------|
| `parallel.py` | Thread/process pool helpers with resource-aware worker sizing |
| `workflow.py` | Config validation, sample config creation, and workflow orchestration |
| `discovery.py` | AST-based function/config discovery, call graphs, and symbol usage search |

## Key Functions and Classes

| Symbol | Description |
|--------|-------------|
| `thread_map()` | Map a function over items using a thread pool with progress |
| `process_map()` | Map a function over items using a process pool |
| `resource_aware_workers()` | Choose worker count based on CPU and memory |
| `rate_limited_map()` | Execute with rate limiting for API-bound workloads |
| `parallel_batch()` | Split items into batches and process in parallel |
| `validate_config_file()` | Validate a config file or dict: a `steps` section, or `downloads` and `processing` sections; optional JSON schema |
| `create_sample_config()` | Write a sample workflow config (`basic`, `scientific` or `advanced`) |
| `download_and_process_data()` | Fetch the files in a config's `downloads` section one at a time. `processing` entries are recorded as `not_executed` and listed in `errors`; the result has no `success` key |
| `run_config_based_workflow()` | Validate and run a config file. A `steps` section runs through `BaseWorkflowOrchestrator`; a `downloads`/`processing` config runs through `download_and_process_data` and the result gains `success`, which is `True` only when `errors` is empty |
| `WorkflowStep` | A named step with its callable, config, dependencies and timing |
| `BaseWorkflowOrchestrator` | Runs `steps` one at a time in dependency order and skips the dependents of a failed step |
| `discover_functions()` | Index every function defined in one module file (private functions and methods included) via AST parsing, with an optional name regex |
| `discover_configs()` | Find config files for a domain across the repository |
| `build_call_graph()` | Trace function calls from an entry point |
| `find_symbol_usage()` | Search for all references to a symbol across the codebase |

## Usage

```python
import json
from pathlib import Path

from metainformant.core.execution.discovery import discover_functions
from metainformant.core.execution.parallel import resource_aware_workers, thread_map
from metainformant.core.execution.workflow import download_and_process_data, run_config_based_workflow

squares = thread_map(lambda n: n * n, range(8), max_workers=resource_aware_workers())

# A `steps` config runs each step's function in dependency order
Path("output/execution_demo").mkdir(parents=True, exist_ok=True)
steps = {"steps": {"cwd": {"function": "os:getcwd"}, "python": {"function": "platform:python_version"}}}
Path("output/execution_demo/steps.json").write_text(json.dumps(steps))
outcome = run_config_based_workflow("output/execution_demo/steps.json")
print(outcome["success"], outcome["execution_order"])  # True ['cwd', 'python']

# A `downloads`/`processing` config fetches the files; processing entries are recorded as not_executed
source = Path("output/execution_demo/source.json")
source.write_text('{"genes": 3}')
results = download_and_process_data(
    {
        "downloads": {"genes": {"url": source.resolve().as_uri(), "filename": "genes.json"}},
        "processing": {"normalization": {"type": "normalize"}},
    },
    output_dir="output/execution_demo/run",
)
print(results["downloads"]["genes"]["status"], results["processing"]["normalization"]["status"])  # success not_executed

# discover_functions takes one module file; `pattern` is an optional regex on function names
funcs = discover_functions("src/metainformant/core/execution/workflow.py", pattern="^download")
print([f.name for f in funcs])  # ['download_and_process_data']
```
