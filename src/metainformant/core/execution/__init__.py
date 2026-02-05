"""Core execution utilities for METAINFORMANT.

This module provides workflow orchestration, parallel processing, and code discovery utilities.
"""

from __future__ import annotations

# Discovery utilities
from .discovery import (
    ConfigInfo,
    FunctionInfo,
    ModuleDependency,
    OutputPattern,
    SymbolUsage,
    build_call_graph,
    discover_configs,
    discover_functions,
    discover_output_patterns,
    discover_workflows,
    find_symbol_usage,
    get_module_dependencies,
    invalidate_discovery_cache,
)

# Parallel processing
from .parallel import (
    cpu_count,
    gather_results,
    parallel_batch,
    process_map,
    rate_limited_map,
    resource_aware_workers,
    thread_map,
    thread_map_unordered,
)

# Workflow orchestration
from .workflow import (
    BaseWorkflowOrchestrator,
    WorkflowStep,
    create_sample_config,
    download_and_process_data,
    run_config_based_workflow,
    validate_config_file,
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
