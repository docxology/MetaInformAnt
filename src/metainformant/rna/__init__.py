"""RNA domain functionality.

This module exposes a thin, modular wrapper around the external `amalgkit`
CLI to support transcriptomic meta-analysis workflows.
"""

from .amalgkit import (
    AmalgkitParams,
    build_amalgkit_command,
    build_cli_args,
    check_cli_available,
    config,
    csca,
    cstmm,
    curate,
    ensure_cli_available,
    getfastq,
    integrate,
    merge,
    metadata,
    quant,
    run_amalgkit,
    sanity,
    select,
)

# Lazy workflow imports to avoid import-time failures on older Python during
# light-weight uses that only need the thin CLI wrappers. Full workflow
# functionality requires Python 3.11+ per project configuration.
try:  # pragma: no cover - exercised in integration tests under py311
    from .workflow import AmalgkitWorkflowConfig, execute_workflow, plan_workflow

    _HAS_WORKFLOW = True
except Exception:  # pragma: no cover - defensive for environments <3.11
    _HAS_WORKFLOW = False

__all__ = [
    "AmalgkitParams",
    "build_cli_args",
    "build_amalgkit_command",
    "check_cli_available",
    "ensure_cli_available",
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
]

if _HAS_WORKFLOW:
    __all__ += [
        "AmalgkitWorkflowConfig",
        "plan_workflow",
        "execute_workflow",
    ]
