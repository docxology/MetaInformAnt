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

# Monitoring and environment functions (always available)
from .environment import (
    check_amalgkit,
    check_dependencies,
    check_kallisto,
    check_metainformant,
    check_rscript,
    check_sra_toolkit,
    check_virtual_env,
    validate_environment,
)
from .monitoring import (
    analyze_species_status,
    check_active_downloads,
    check_workflow_progress,
    count_quantified_samples,
    find_unquantified_samples,
    get_sample_status,
)

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
    # Monitoring
    "analyze_species_status",
    "check_active_downloads",
    "check_workflow_progress",
    "count_quantified_samples",
    "find_unquantified_samples",
    "get_sample_status",
    # Environment
    "check_amalgkit",
    "check_dependencies",
    "check_kallisto",
    "check_metainformant",
    "check_rscript",
    "check_sra_toolkit",
    "check_virtual_env",
    "validate_environment",
]

if _HAS_WORKFLOW:
    __all__ += [
        "AmalgkitWorkflowConfig",
        "plan_workflow",
        "execute_workflow",
    ]
