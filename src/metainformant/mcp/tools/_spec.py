"""Shared contract for MetaInformAnt MCP-adjacent analysis tools.

Every tool in ``metainformant.mcp.tools`` exposes a module-level
``TOOL_SPEC`` dictionary describing the tool's name, description, JSON-schema
input, and handler entry point. The shape is intentionally minimal so that a
later MCP registry (e.g. the mcp1 lane's ``registry.py``) can enumerate and
register tools by importing ``TOOL_SPEC`` from each module.

Handler contract:
- Accept the documented input fields as keyword arguments.
- Return a JSON-serializable ``dict``.
- Be deterministic for identical inputs.
- Read-only tools touch only paths given as inputs; writing tools create
  files only under a caller-supplied output directory.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

Handler = Callable[..., dict]


def validate_output_dir(path: str | Path) -> Path:
    """Resolve and create an explicit output directory, returning its Path.

    Every writing tool must funnel its output location through this helper so
    the 'explicit-output-dir writes only' invariant holds uniformly.
    """
    resolved = Path(path).expanduser().resolve()
    resolved.mkdir(parents=True, exist_ok=True)
    if not resolved.is_dir():
        raise ValueError(f"output_dir did not resolve to a directory: {resolved}")
    return resolved


def dump_json(obj: Any, path: Path) -> Path:
    """Write deterministic JSON (sorted keys, fixed indent) and return path."""
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def read_table(path: str | Path, sep: str | None = None):
    """Read a delimited expression/summary table into a DataFrame.

    Auto-detects comma vs tab when ``sep`` is None. Raises FileNotFoundError
    (not a silent empty frame) so callers see real input errors.
    """
    import pandas as pd

    p = Path(path).expanduser()
    if not p.exists():
        raise FileNotFoundError(f"input table not found: {p}")
    if sep is None:
        sep = "\t" if p.suffix in {".tsv", ".tab"} else ","
    return pd.read_csv(p, sep=sep, index_col=0)


__all__ = ["Handler", "validate_output_dir", "dump_json", "read_table"]
