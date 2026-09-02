"""Core MCP tools: path resolution and config-file validation (read-only)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from metainformant.core.io.paths import is_safe_path, is_within


def _handle_path_resolve(path: str, parent: str | None = None) -> dict:
    """Resolve a path and optionally check containment within a parent (read-only)."""
    resolved = Path(path).expanduser().resolve()
    result: dict[str, Any] = {
        "path": str(resolved),
        "exists": resolved.exists(),
        "is_file": resolved.is_file(),
        "is_dir": resolved.is_dir(),
        "safe": is_safe_path(path),
    }
    if parent is not None:
        parent_resolved = Path(parent).expanduser().resolve()
        result["within_parent"] = is_within(resolved, parent_resolved)
        result["parent"] = str(parent_resolved)
    return result


def _handle_config_validate(config_path: str) -> dict:
    """Load and validate a JSON/YAML config file (read-only)."""
    p = Path(config_path).expanduser()
    if not p.exists():
        return {"valid": False, "error": f"config not found: {p}"}
    from metainformant.core.utils.config import load_mapping_from_file

    try:
        data = load_mapping_from_file(p)
    except Exception as exc:  # real parse errors reported, not swallowed
        return {"valid": False, "error": f"{type(exc).__name__}: {exc}"}
    if not isinstance(data, dict):
        return {"valid": False, "error": f"top-level structure is {type(data).__name__}, expected mapping"}
    return {"valid": True, "top_level_keys": sorted(data.keys()), "n_keys": len(data)}


TOOL_SPEC: dict[str, Any] = {
    "name": "core_path_resolve",
    "description": "Resolve a filesystem path and check safety/containment (read-only).",
    "input_schema": {
        "type": "object",
        "properties": {"path": {"type": "string"}, "parent": {"type": "string"}},
        "required": ["path"],
    },
    "handler": _handle_path_resolve,
    "writes": "read-only",
}

CONFIG_SPEC: dict[str, Any] = {
    "name": "core_config_validate",
    "description": "Load and validate a JSON/YAML config mapping file (read-only).",
    "input_schema": {
        "type": "object",
        "properties": {"config_path": {"type": "string"}},
        "required": ["config_path"],
    },
    "handler": _handle_config_validate,
    "writes": "read-only",
}

ALL_SPECS: list[dict[str, Any]] = [TOOL_SPEC, CONFIG_SPEC]
