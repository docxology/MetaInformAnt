# Amalgkit Centralized Scripts

This directory contains centralized scripts and utilities for amalgkit workflows that work across **all species** without needing to be copied.

## Scripts

### `run_amalgkit.sh`
**Purpose**: Comprehensive orchestrator for end-to-end amalgkit pipeline  
**Location**: `scripts/rna/amalgkit/run_amalgkit.sh`

**Features**:
- Automatic environment setup (uv, venv, dependencies)
- Intelligent step execution with dependency checking
- Robust FASTQ downloading with multiple fallback sources
- Comprehensive logging and monitoring
- Progress tracking and resource monitoring

**Usage Options**:
```bash
# Full pipeline run
scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit_pogonomyrmex_barbatus.yaml

# Run specific steps only
scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit_pogonomyrmex_barbatus.yaml --steps metadata,select,getfastq

# Skip expensive steps for testing
scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit_pogonomyrmex_barbatus.yaml --skip-steps getfastq,quant

# Dry run to see execution plan
scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit_pogonomyrmex_barbatus.yaml --dry-run
```

### `verify_workflow.sh`
**Purpose**: Validates amalgkit workflow outputs for any species  
**Location**: `scripts/rna/amalgkit/verify_workflow.sh`

**Usage Options**:

```bash
# Option 1: From repo root with species name
scripts/rna/amalgkit/verify_workflow.sh pbarbatus
scripts/rna/amalgkit/verify_workflow.sh sinvicta

# Option 2: From species directory
cd output/amalgkit/pbarbatus
../../../scripts/rna/amalgkit/verify_workflow.sh .

# Option 3: With full path
scripts/rna/amalgkit/verify_workflow.sh output/amalgkit/cfloridanus
```

**What it checks**:
- ✅ Data integrity with `amalgkit sanity`
- ✅ Curate outputs (6 PDFs, 7 TSVs, 3 RData)
- ✅ Sample counts and file completeness
- ✅ Auto-regenerates curate if incomplete
- ✅ Provides workflow status summary

**Output**: Comprehensive verification report with file counts and quick access commands

### `run_workflow.py` (Main Orchestrator)
**Purpose**: Orchestrates full workflows for single species (run separately for multiple species)  
**Location**: `scripts/rna/run_workflow.py`

**Features**:
- Complete end-to-end workflow execution
- All 11 amalgkit steps (metadata → sanity)
- Status checking and cleanup operations
- Parallel downloads via `num_download_workers` configuration

**Usage**:
```bash
# Single species
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Multiple species (run separately for each)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species1.yaml
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_species2.yaml
```

### Future Scripts
- `compare_species.sh` - Compare expression across multiple species
- `extract_top_genes.sh` - Extract highly expressed genes  
- `validate_cross_species.sh` - Validate cross-species analysis outputs

## Why Centralized?


### After (Clean)
```
scripts/rna/amalgkit/
└── verify_workflow.sh                 ✅ Single source of truth

output/amalgkit/
├── pbarbatus/work/                    ✅ Only data outputs
├── sinvicta/work/                     ✅ Only data outputs
└── cfloridanus/work/                  ✅ Only data outputs
```
- ✅ Single script works for all species
- ✅ Tracked in git, versioned properly
- ✅ Update once, applies everywhere
- ✅ Output directories stay clean (data only)

## Integration

These scripts are referenced in:
- `docs/rna/GETTING_STARTED.md` - Complete setup and workflow guides
- `docs/rna/amalgkit/README.md` - Overview and links
- `output/amalgkit/COMPREHENSIVE_WORKFLOW_STATUS.md` - Execution status (program-generated)

## Usage Pattern

1. **Run script directly** from scripts/ directory
2. **Point to any species** by name or path
3. **Review results** and follow recommended next steps
4. **No copying needed** - scripts stay in one place

---

**Location**: `scripts/rna/amalgkit/`  
**Documentation**: `docs/rna/amalgkit/`  
**Philosophy**: Centralized methods, distributed data

