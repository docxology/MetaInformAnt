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
| `validate_config_file()` | Validate a YAML config against expected schema |
| `run_config_based_workflow()` | Execute a full workflow from a config file |
| `WorkflowStep` | Dataclass defining a named step with its callable |
| `BaseWorkflowOrchestrator` | Base class for multi-step workflow execution |
| `discover_functions()` | Index all public functions in a module via AST parsing |
| `discover_configs()` | Find config files for a domain across the repository |
| `build_call_graph()` | Trace function calls from an entry point |
| `find_symbol_usage()` | Search for all references to a symbol across the codebase |

## Usage

```python
from metainformant.core.execution.parallel import thread_map, resource_aware_workers
from metainformant.core.execution.workflow import run_config_based_workflow
from metainformant.core.execution.discovery import discover_functions

results = thread_map(process_sample, samples, max_workers=resource_aware_workers())
run_config_based_workflow("config/amalgkit/workflow.yaml")
funcs = discover_functions("src/metainformant/dna/")
```
