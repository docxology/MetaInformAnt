"""RNA workflow orchestration engine.

This subpackage provides the workflow execution engine including:
- Workflow planning, execution, and validation
- Species discovery via NCBI
- Progress tracking and monitoring
- Pipeline summarization
- High-level orchestration
"""

from __future__ import annotations

from .discovery import (
    generate_config_yaml,
    get_genome_info,
    search_species_with_rnaseq,
)
from .monitoring import (
    analyze_species_status,
    check_workflow_progress,
    count_quantified_samples,
    find_unquantified_samples,
    get_sample_progress_report,
)
from .orchestration import (
    cleanup_unquantified_samples,
    discover_species_configs,
    estimate_workflow_resources,
    get_workflow_status,
    monitor_workflows,
    retry_failed_steps,
    run_parallel_workflows,
    run_workflow_for_species,
    validate_and_execute,
)
from .pipeline import (
    summarize_curate_tables,
)
from .progress_tracker import (
    ProgressTracker,
    get_tracker,
)
from .sra_extraction import (
    extract_sra_directly,
    manual_integration_fallback,
)
from .workflow import (
    AmalgkitWorkflowConfig,
    WorkflowExecutionResult,
    WorkflowStepResult,
    create_sample_config,
    execute_workflow,
    load_workflow_config,
    plan_workflow,
    run_config_based_workflow,
    validate_workflow_config,
    validate_workflow_outputs,
)
from .workflow_cleanup import (
    check_disk_space,
    check_disk_space_or_fail,
    cleanup_after_quant,
    cleanup_fastqs,
    cleanup_incorrectly_placed_sra_files,
    cleanup_temp_files,
    filter_metadata_for_unquantified,
    get_quantified_samples,
)
from .workflow_steps import (
    check_step_completion_status,
    handle_post_step_actions,
    log_workflow_summary,
    setup_vdb_config,
    validate_step_prerequisites,
)

__all__ = [
    # Workflow core
    "AmalgkitWorkflowConfig",
    "WorkflowStepResult",
    "WorkflowExecutionResult",
    "load_workflow_config",
    "execute_workflow",
    "plan_workflow",
    "validate_workflow_config",
    "validate_workflow_outputs",
    "run_config_based_workflow",
    "create_sample_config",
    # Progress
    "ProgressTracker",
    "get_tracker",
    # Pipeline
    "summarize_curate_tables",
    # Orchestration
    "run_workflow_for_species",
    "cleanup_unquantified_samples",
    "monitor_workflows",
    "discover_species_configs",
    "run_parallel_workflows",
    "validate_and_execute",
    "retry_failed_steps",
    "get_workflow_status",
    "estimate_workflow_resources",
    # Monitoring
    "check_workflow_progress",
    "analyze_species_status",
    "count_quantified_samples",
    "find_unquantified_samples",
    "get_sample_progress_report",
    # Discovery
    "search_species_with_rnaseq",
    "get_genome_info",
    "generate_config_yaml",
    # Workflow cleanup
    "check_disk_space",
    "check_disk_space_or_fail",
    "cleanup_temp_files",
    "cleanup_incorrectly_placed_sra_files",
    "cleanup_fastqs",
    "get_quantified_samples",
    "cleanup_after_quant",
    "filter_metadata_for_unquantified",
    # SRA extraction
    "extract_sra_directly",
    "manual_integration_fallback",
    # Workflow steps
    "setup_vdb_config",
    "check_step_completion_status",
    "validate_step_prerequisites",
    "handle_post_step_actions",
    "log_workflow_summary",
]
