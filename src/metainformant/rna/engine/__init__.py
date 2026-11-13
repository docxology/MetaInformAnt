"""Current RNA/Amalgkit execution engine.

The production path is the streaming producer backed by SQLite progress and
the hash-bound downstream checkpoint runner. The workflow modules remain
internal planning and command-construction helpers used by that producer.

Submodules load on first attribute access (PEP 562). Importing the package
therefore does not import ``streaming_orchestrator`` and, through it,
``campaign_lock`` and ``preflight``, so ``python -m
metainformant.rna.engine.<module>`` runs without runpy's "found in
sys.modules" warning.
"""

from __future__ import annotations

import importlib
from types import ModuleType

_SUBMODULES = (
    "discovery",
    "evidence_bundle",
    "pipeline",
    "process_observation",
    "progress_db",
    "sra_extraction",
    "streaming_orchestrator",
    "workflow",
    "workflow_cleanup",
    "workflow_core",
    "workflow_execution",
    "workflow_planning",
    "workflow_steps",
)

__all__ = list(_SUBMODULES)


def __getattr__(name: str) -> ModuleType:
    if name in _SUBMODULES:
        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(_SUBMODULES))
