"""RNA-seq workflow orchestration and configuration management.

This module provides high-level functions for managing RNA-seq analysis workflows,
including configuration loading, workflow planning, and execution orchestration.

This is a re-export module that aggregates the public API from:
- workflow_core: Config classes, validation, sample config creation
- workflow_planning: Workflow planning, step defaults, reference genome prep
- workflow_execution: Workflow execution engine, streaming mode
"""

from __future__ import annotations

# Re-export core classes and config functions
from metainformant.rna.engine.workflow_core import (
    AmalgkitWorkflowConfig,
    WorkflowExecutionResult,
    WorkflowStepResult,
    apply_config_defaults,
    create_sample_config,
    load_workflow_config,
    validate_workflow_config,
    validate_workflow_outputs,
)

# Re-export planning functions
from metainformant.rna.engine.workflow_planning import (
    _is_step_completed,
    _log_getfastq_summary,
    _log_heartbeat,
    apply_step_defaults,
    create_extraction_metadata,
    plan_workflow,
    plan_workflow_with_params,
    prepare_extraction_directories,
    prepare_reference_genome,
    sanitize_params_for_cli,
    verify_getfastq_prerequisites,
)

# Re-export execution functions
from metainformant.rna.engine.workflow_execution import (
    execute_workflow,
    run_config_based_workflow,
)

# Backward compatibility aliases
from metainformant.rna.engine.workflow_cleanup import (
    check_disk_space,
    check_disk_space_or_fail,
    cleanup_incorrectly_placed_sra_files,
    cleanup_temp_files,
)

_cleanup_incorrectly_placed_sra_files = cleanup_incorrectly_placed_sra_files
_cleanup_temp_files = cleanup_temp_files
_check_disk_space = check_disk_space
_check_disk_space_or_fail = check_disk_space_or_fail

__all__ = [
    # Core classes
    "WorkflowStepResult",
    "WorkflowExecutionResult",
    "AmalgkitWorkflowConfig",
    # Config functions
    "load_workflow_config",
    "apply_config_defaults",
    "apply_step_defaults",
    "validate_workflow_config",
    "validate_workflow_outputs",
    "create_sample_config",
    # Planning functions
    "plan_workflow",
    "plan_workflow_with_params",
    "prepare_reference_genome",
    "prepare_extraction_directories",
    "create_extraction_metadata",
    "verify_getfastq_prerequisites",
    "sanitize_params_for_cli",
    # Execution functions
    "execute_workflow",
    "run_config_based_workflow",
]
