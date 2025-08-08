"""RNA domain functionality.

This module exposes a thin, modular wrapper around the external `amalgkit`
CLI to support transcriptomic meta-analysis workflows.
"""

from .amalgkit import (
    AmalgkitParams,
    build_amalgkit_command,
    build_cli_args,
    check_cli_available,
    run_amalgkit,
    metadata,
    integrate,
    config,
    select,
    getfastq,
    quant,
    merge,
    cstmm,
    curate,
    csca,
    sanity,
)
from .workflow import (
    AmalgkitWorkflowConfig,
    plan_workflow,
    execute_workflow,
)

__all__ = [
    "AmalgkitParams",
    "build_cli_args",
    "build_amalgkit_command",
    "check_cli_available",
    "run_amalgkit",
    "metadata",
    "integrate",
    "config",
    "select",
    "getfastq",
    "quant",
    "merge",
    "cstmm",
    "curate",
    "csca",
    "sanity",
    "AmalgkitWorkflowConfig",
    "plan_workflow",
    "execute_workflow",
]



