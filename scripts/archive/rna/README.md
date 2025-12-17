# Archived RNA Scripts

This directory contains shell scripts and Python scripts that have been archived because their functionality is now provided by `run_workflow.py` or Python modules.

## Archived Scripts

### check_getfastq_status.sh
**Status**: Archived  
**Reason**: Functionality provided by `run_workflow.py --status --detailed` and `metainformant.rna.monitoring.check_workflow_progress()`  
**Replacement**: `python3 scripts/rna/run_workflow.py --config <config> --status --detailed`

### monitor_getfastq.sh
**Status**: Archived  
**Reason**: Wrapper around `check_getfastq_status.sh` (also archived)  
**Replacement**: Use `run_workflow.py --status` in a watch loop, or use Python monitoring functions

### monitor_progress.sh
**Status**: Archived  
**Reason**: Wrapper around `check_getfastq_status.sh` (also archived)  
**Replacement**: Use `run_workflow.py --status` in a watch loop, or use Python monitoring functions

### view_dashboard.sh
**Status**: Archived  
**Reason**: Checks for dashboard files that are no longer generated  
**Replacement**: Use `run_workflow.py --status` for status information

### quant_and_cleanup.sh
**Status**: Archived  
**Reason**: Functionality provided by `run_workflow.py --cleanup-unquantified` and `metainformant.rna.orchestration.cleanup_unquantified_samples()`  
**Replacement**: `python3 scripts/rna/run_workflow.py --config <config> --cleanup-unquantified`

### run_quant_cleanup.py
**Status**: Archived  
**Reason**: Functionality provided by `run_workflow.py --cleanup-unquantified`  
**Replacement**: `python3 scripts/rna/run_workflow.py --config <config> --cleanup-unquantified`

### run_pogonomyrmex_barbatus_workflow.sh
**Status**: Archived  
**Reason**: Wrapper script that calls `run_workflow.py` multiple times. Can be replaced with direct `run_workflow.py` calls  
**Replacement**: `python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml`

## Still Active Scripts

The following scripts remain in `scripts/rna/` because they provide unique functionality:

- **fix_tmp_space.sh**: Utility for cleaning up /tmp space (moved to `scripts/package/`)
- **setup_genome.py**: Genome setup utility (unique functionality)
- **discover_species.py**: Species discovery utility (unique functionality)
- **check_environment.py**: Environment checking utility (unique functionality)
- **run_workflow.py**: Main workflow orchestrator (recommended for all workflows)

