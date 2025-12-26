# RNA Module Architecture

This document provides a comprehensive overview of the METAINFORMANT RNA module architecture, explaining how components interact and how to extend the system.

## Table of Contents

1. [Module Structure](#module-structure)
2. [Processing Workflows](#processing-workflows)
3. [Step Runner Pattern](#step-runner-pattern)
4. [Integration Points](#integration-points)
5. [Data Flow](#data-flow)
6. [Decision Tree](#decision-tree)

## Module Structure

### Directory Layout

```
src/metainformant/rna/
├── __init__.py              # Public API exports
├── workflow.py              # Main workflow orchestration
├── orchestration.py         # High-level workflow execution
├── amalgkit.py              # Amalgkit CLI wrapper
├── genome_prep.py           # Genome preparation utilities
├── configs.py               # Configuration management
├── steps/                   # Step implementations
│   ├── __init__.py          # Step runner registry
│   ├── metadata.py          # Step 1: Metadata retrieval
│   ├── config.py            # Step 2: Config generation
│   ├── select.py            # Step 3: Sample selection
│   ├── getfastq.py          # Step 4: FASTQ download
│   ├── integrate.py         # Step 5: FASTQ integration
│   ├── quant.py             # Step 6: Quantification
│   ├── merge.py             # Step 7: Result merging
│   ├── cstmm.py             # Step 8: Cross-species normalization
│   ├── curate.py            # Step 9: Quality curation
│   ├── csca.py              # Step 10: Correlation analysis
│   ├── sanity.py            # Step 11: Final validation
│   ├── process_samples.py   # Unified download-quantify workflow
│   └── download_progress.py # Progress tracking
└── ...
```

### Component Responsibilities

#### Core Orchestration
- **`workflow.py`**: Plans and executes complete workflows based on configuration
- **`orchestration.py`**: High-level functions for running workflows for species
- **`amalgkit.py`**: Low-level wrapper for `amalgkit` CLI commands

#### Step Runners
- **11 step modules**: Thin wrappers around `amalgkit` subcommands
- **`process_samples.py`**: Unified workflow for download-quantify-delete operations
- Each step is independent and can be run standalone or as part of a workflow

#### Utilities
- **`genome_prep.py`**: Downloads and indexes reference genomes
- **`configs.py`**: Loads and validates workflow configurations
- **`download_progress.py`**: Real-time progress tracking

## Processing Workflows

The RNA module supports two processing modes for handling large cohorts:

### Sequential Mode (`num_workers=1`)

**Architecture**: One sample at a time
```
For each sample:
  1. Download FASTQ (getfastq)
  2. Wait for conversion to complete
  3. Quantify (quant)
  4. Delete FASTQ files
  5. Next sample
```

**Characteristics**:
- Maximum disk efficiency: only one sample's FASTQs exist at a time
- Lower throughput: samples processed sequentially
- Best for: Limited disk space, small cohorts, or when disk space is critical

**Use Case**: When disk space is very limited or you want maximum control over resource usage.

### Parallel Mode (`num_workers>1`)

**Architecture**: Producer-consumer pattern
```
N download workers (parallel):
  - Download FASTQ files concurrently
  - Place completed downloads in queue

1 quantification worker (sequential):
  - Processes downloads from queue
  - Quantifies one at a time
  - Deletes FASTQ after quantification
```

**Characteristics**:
- Higher throughput: multiple downloads in parallel
- Sequential quantification prevents disk exhaustion
- Best for: Large cohorts, sufficient disk space, faster processing

**Use Case**: When you have sufficient disk space and want to maximize download throughput.

### Unified Implementation

Both modes are implemented in a single function: `run_download_quant_workflow()` in `process_samples.py`. The `num_workers` parameter controls the mode:

```python
from metainformant.rna.steps import run_download_quant_workflow

# Sequential mode
stats = run_download_quant_workflow(
    metadata_path="metadata.tsv",
    getfastq_params={...},
    quant_params={...},
    num_workers=1,  # Sequential
)

# Parallel mode
stats = run_download_quant_workflow(
    metadata_path="metadata.tsv",
    getfastq_params={...},
    quant_params={...},
    num_workers=8,  # 8 parallel download workers
)
```

## Step Runner Pattern

### Design Principles

1. **Thin Wrappers**: Step modules are thin wrappers around `amalgkit` CLI commands
2. **Uniform Interface**: All steps follow the same function signature
3. **Isolation**: Step-specific logic is isolated from orchestration
4. **Composability**: Steps can be run individually or as part of workflows

### Step Interface

All step runners follow this pattern:

```python
def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run amalgkit step.
    
    Args:
        params: Step-specific parameters
        work_dir: Working directory
        log_dir: Log directory
        check: Raise on failure
        
    Returns:
        CompletedProcess with return code and output
    """
    from ..amalgkit import step_name as _step
    return _step(params, work_dir=work_dir, log_dir=log_dir, step_name="step_name", check=check)
```

### Step Registry

Steps are registered in `steps/__init__.py`:

```python
STEP_RUNNERS: dict[str, Callable[..., object]] = {
    "metadata": run_metadata,
    "config": run_config,
    "select": run_select,
    "getfastq": run_getfastq,
    "integrate": run_integrate,
    "quant": run_quant,
    "merge": run_merge,
    "cstmm": run_cstmm,
    "curate": run_curate,
    "csca": run_csca,
    "sanity": run_sanity,
}
```

This registry allows workflows to dynamically call steps by name.

## Integration Points

### Scripts → Source → Amalgkit

```
scripts/rna/run_workflow.py
    ↓
src/metainformant/rna/workflow.py
    ↓
src/metainformant/rna/steps/*.py
    ↓
src/metainformant/rna/amalgkit.py
    ↓
amalgkit CLI (subprocess)
```

### Configuration Flow

```
config/amalgkit/species.yaml
    ↓
metainformant.core.config.load_workflow_config()
    ↓
AmalgkitWorkflowConfig
    ↓
workflow.py:plan_workflow() / execute_workflow()
    ↓
Step runners
```

### Data Flow

```
Metadata TSV
    ↓
Step 1: metadata → metadata.tsv
    ↓
Step 2: config → .config files
    ↓
Step 3: select → filtered metadata
    ↓
Step 4: getfastq → FASTQ files
    ↓
Step 5: integrate → metadata with FASTQ paths
    ↓
Step 6: quant → abundance.tsv per sample
    ↓
Step 7: merge → expression matrix
    ↓
Step 8-11: Statistical analysis
```

## Data Flow

### Workflow Execution Flow

```
1. Load Configuration
   └─> config/amalgkit/species.yaml
       └─> AmalgkitWorkflowConfig

2. Plan Workflow
   └─> workflow.py:plan_workflow()
       └─> Returns list of step names

3. Execute Steps
   └─> workflow.py:execute_workflow()
       ├─> For each step:
       │   ├─> Get step runner from STEP_RUNNERS
       │   ├─> Prepare step parameters
       │   ├─> Call step runner
       │   └─> Handle errors
       │
       └─> Special handling for getfastq+quant:
           └─> process_samples.py:run_download_quant_workflow()
               ├─> Sequential mode (num_workers=1)
               └─> Parallel mode (num_workers>1)
```

### File System Structure

```
output/amalgkit/<species>/
├── work/
│   ├── metadata/
│   │   └── metadata.tsv
│   ├── config/
│   │   └── *.config
│   └── index/
│       └── kallisto.idx
├── fastq/
│   └── getfastq/
│       └── <SRR>/
│           └── *.fastq.gz
├── quant/
│   └── <SRR>/
│       └── abundance.tsv
├── merged/
│   └── expression_matrix.tsv
└── logs/
    └── *.log
```

## Decision Tree

### When to Use Sequential Mode

```
Do you have limited disk space?
├─> Yes → Use sequential mode (num_workers=1)
│   └─> Only one sample's FASTQs exist at a time
│
└─> No → Continue to next question

Is your cohort small (<50 samples)?
├─> Yes → Sequential mode is fine
│   └─> Simpler, easier to debug
│
└─> No → Consider parallel mode

Do you need maximum throughput?
├─> Yes → Use parallel mode (num_workers=4-8)
│   └─> Multiple downloads in parallel
│
└─> No → Sequential mode is sufficient
```

### When to Use Parallel Mode

```
Do you have sufficient disk space?
├─> No → Use sequential mode
│
└─> Yes → Continue

Is your cohort large (>50 samples)?
├─> No → Sequential mode is fine
│
└─> Yes → Use parallel mode

How many workers?
├─> Limited bandwidth → 2-4 workers
├─> Good bandwidth → 4-8 workers
└─> Excellent bandwidth → 8-16 workers
    (Note: More workers = more disk space needed)
```

### Step Selection

```
What do you want to do?
├─> Start from scratch
│   └─> Steps: metadata, config, select, getfastq, quant, merge
│
├─> Already have metadata
│   └─> Steps: config, select, getfastq, quant, merge
│
├─> Already have FASTQ files
│   └─> Steps: integrate, quant, merge
│
├─> Already have quantifications
│   └─> Steps: merge, cstmm, curate, csca, sanity
│
└─> Complete workflow
    └─> Steps: metadata, config, select, getfastq, quant, merge, cstmm, curate, csca, sanity
```

## Extending the System

### Adding a New Step

1. **Create step module** (`steps/new_step.py`):
```python
from ..amalgkit import new_step as _new_step

def run(params, *, work_dir=None, log_dir=None, check=False):
    """Run amalgkit new_step."""
    return _new_step(params, work_dir=work_dir, log_dir=log_dir, step_name="new_step", check=check)
```

2. **Register in `steps/__init__.py`**:
```python
from .new_step import run as run_new_step

STEP_RUNNERS["new_step"] = run_new_step
__all__.append("run_new_step")
```

3. **Add to workflow planning** (if needed):
```python
# In workflow.py:plan_workflow()
if "new_step" in config.steps:
    planned_steps.append("new_step")
```

4. **Add documentation**:
   - Update `docs/rna/amalgkit/steps/README.md`
   - Create `docs/rna/amalgkit/steps/new_step.md`

### Modifying Processing Logic

The unified processing function (`process_samples.py`) can be extended:

- **Custom progress tracking**: Pass custom `DownloadProgressMonitor`
- **Custom error handling**: Modify worker functions
- **Custom cleanup**: Modify `_delete_fastq_for_sample()`

## Related Documentation

- [Workflow Guide](workflow.md) - Workflow planning and execution
- [Step Documentation](amalgkit/steps/README.md) - Individual step guides
- [Configuration Guide](CONFIGURATION.md) - Configuration management
- [API Reference](API.md) - Complete function documentation

