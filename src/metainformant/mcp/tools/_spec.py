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
import math
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd

Handler = Callable[..., dict[str, Any]]


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


def json_safe_float(value: Any) -> float | None:
    """Return ``value`` as a plain float, with NaN/inf mapped to None.

    Strict JSON has no NaN/Infinity tokens, so descriptive statistics over
    degenerate inputs (e.g. std of a single row) must surface as null.
    """
    number = float(value)
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def dump_json(obj: Any, path: Path) -> Path:
    """Write deterministic JSON (sorted keys, fixed indent) and return path.

    ``allow_nan=False`` makes any non-JSON-safe float a loud error instead of
    a silently invalid artifact; use :func:`json_safe_float` when converting
    possibly-degenerate statistics.
    """
    payload = json.dumps(obj, indent=2, sort_keys=True, allow_nan=False)
    path.write_text(payload + "\n", encoding="utf-8")
    return path


def read_table(
    path: str | Path, sep: str | None = None, index_col: int | str | Sequence[int | str] | None = 0
) -> pd.DataFrame:
    """Read a delimited expression/summary table into a DataFrame.
    Auto-detects comma vs tab when ``sep`` is None. Pass index_col=None for
    tables whose first column is plain data. Raises FileNotFoundError
    (not a silent empty frame) so callers see real input errors.
    """
    import pandas as pd

    p = Path(path).expanduser()
    if not p.exists():
        raise FileNotFoundError(f"input table not found: {p}")
    if sep is None:
        sep = "\t" if p.suffix in {".tsv", ".tab"} else ","
    return pd.read_csv(p, sep=sep, index_col=index_col)


__all__ = ["Handler", "validate_output_dir", "dump_json", "json_safe_float", "read_table"]
