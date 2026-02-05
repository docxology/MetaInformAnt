"""Core execution utilities for METAINFORMANT.

This module provides workflow orchestration, parallel processing, and code discovery utilities.
"""

from __future__ import annotations

# Discovery utilities
from .discovery import (
    FunctionInfo,
    ConfigInfo,
    OutputPattern,
    SymbolUsage,
    ModuleDependency,
    discover_functions,
    discover_configs,
    discover_output_patterns,
    build_call_graph,
    find_symbol_usage,
    get_module_dependencies,
    discover_workflows,
    invalidate_discovery_cache,
)

# Parallel processing
from .parallel import (
    thread_map,
    thread_map_unordered,
    process_map,
    parallel_batch,
    cpu_count,
    gather_results,
    rate_limited_map,
    resource_aware_workers,
)

# Workflow orchestration
from .workflow import (
    validate_config_file,
    create_sample_config,
    download_and_process_data,
    run_config_based_workflow,
    WorkflowStep,
    BaseWorkflowOrchestrator,
)

__all__ = [
    # Discovery
    "FunctionInfo",
    "ConfigInfo",
    "OutputPattern",
    "SymbolUsage",
    "ModuleDependency",
    "discover_functions",
    "discover_configs",
    "discover_output_patterns",
    "build_call_graph",
    "find_symbol_usage",
    "get_module_dependencies",
    "discover_workflows",
    "invalidate_discovery_cache",
    # Parallel
    "thread_map",
    "thread_map_unordered",
    "process_map",
    "parallel_batch",
    "cpu_count",
    "gather_results",
    "rate_limited_map",
    "resource_aware_workers",
    # Workflow
    "validate_config_file",
    "create_sample_config",
    "download_and_process_data",
    "run_config_based_workflow",
    "WorkflowStep",
    "BaseWorkflowOrchestrator",
]
