"""Tool adapter wiring the existing amalgkit monitor into the MCP registry.

This module exposes the pre-existing
:mod:`metainformant.mcp.tools.amalgkit_monitor` status snapshot as a declarative
MCP tool.  The monitor module itself is imported, never modified.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Mapping

from metainformant.mcp.registry import Tool
from metainformant.mcp.tools import amalgkit_monitor

_DEFAULT_DATA_ROOT = Path(os.environ.get("AMALGKIT_DATA_ROOT", "output/amalgkit"))


def _run_amalgkit_monitor(arguments: Mapping[str, Any]) -> dict[str, Any]:
    """Execute the monitor snapshot with validated, schema-typed arguments."""

    data_root = Path(str(arguments.get("data_root") or _DEFAULT_DATA_ROOT))
    log_file_raw = arguments.get("log_file")
    log_file = Path(str(log_file_raw)) if log_file_raw else None
    inspect_processes = bool(arguments.get("inspect_processes", False))
    return amalgkit_monitor.build_status(
        data_root=data_root,
        log_file=log_file,
        inspect_processes=inspect_processes,
    )


TOOL = Tool(
    name="amalgkit_monitor",
    description=(
        "Machine-readable operational snapshot of the Amalgkit RNA campaign: "
        "process status, log-derived progress, database-backed cohort readiness "
        "(descriptive only; biological inference always withheld), and system stats."
    ),
    input_schema={
        "type": "object",
        "properties": {
            "data_root": {
                "type": "string",
                "description": (
                    "Campaign data root containing pipeline_progress.db " f"(default: {_DEFAULT_DATA_ROOT})."
                ),
            },
            "log_file": {
                "type": "string",
                "description": "Optional explicit path to the campaign log file.",
            },
            "inspect_processes": {
                "type": "boolean",
                "description": (
                    "Scan the process table for running workflow processes "
                    "(default false; locks and receipts are used instead)."
                ),
            },
        },
        "required": [],
    },
    handler=_run_amalgkit_monitor,
    metadata={"module": "metainformant.mcp.tools.amalgkit_monitor", "read_only": True},
)

TOOLS = [TOOL]
