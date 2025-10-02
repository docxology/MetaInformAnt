from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Mapping

# Lightweight, optional YAML/TOML support without hard deps
try:  # pragma: no cover - optional dependency
    import yaml  # type: ignore
except Exception:  # pragma: no cover - optional
    yaml = None  # type: ignore

try:  # pragma: no cover - optional dependency
    import importlib

    tomllib = importlib.import_module("tomllib")  # Python 3.11+
except Exception:  # pragma: no cover - optional
    tomllib = None  # type: ignore


@dataclass(frozen=True)
class PostgresConfig:
    host: str
    port: int
    database: str
    user: str
    password: str


def load_postgres_config_from_env(prefix: str = "PG") -> PostgresConfig | None:
    host = os.getenv(f"{prefix}_HOST")
    database = os.getenv("DB_NAME", os.getenv(f"{prefix}_DATABASE"))
    user = os.getenv("DB_USER", os.getenv(f"{prefix}_USER"))
    password = os.getenv("DB_PASSWORD", os.getenv(f"{prefix}_PASSWORD"))
    port_str = os.getenv(f"{prefix}_PORT", "5432")

    if not host or not database or not user or not password:
        return None

    try:
        port = int(port_str)
    except ValueError:
        port = 5432

    return PostgresConfig(host=host, port=port, database=database, user=user, password=password)


def load_config_file(config_path: Path) -> Dict[str, Any]:
    """Load configuration from YAML, TOML, or JSON file.

    Args:
        config_path: Path to configuration file

    Returns:
        Configuration dictionary

    Raises:
        ValueError: If file format not supported or file cannot be parsed
    """
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    suffix = config_path.suffix.lower()
    content = config_path.read_text()

    if suffix in {".yaml", ".yml"}:
        if yaml is None:
            raise RuntimeError("PyYAML not available for YAML config files")
        return yaml.safe_load(content)
    elif suffix == ".toml":
        if tomllib is None:
            raise RuntimeError("tomllib not available for TOML config files")
        return tomllib.loads(content)
    elif suffix == ".json":
        import json

        return json.loads(content)
    else:
        raise ValueError(f"Unsupported config file format: {suffix}")


def get_env_or_default(env_var: str, default: str) -> str:
    """Get environment variable value or default.

    Args:
        env_var: Environment variable name
        default: Default value if not set

    Returns:
        Environment variable value or default
    """
    return os.getenv(env_var, default)


def _coerce_bool(s: str) -> bool:
    return s.strip().lower() in {"1", "true", "yes", "y", "on"}


def load_typed_env(*, prefix: str, keys: Mapping[str, type]) -> Dict[str, Any]:
    """Load a mapping of env vars with types.

    Example:
        load_typed_env(prefix="APP", keys={"PORT": int, "DEBUG": bool})
    """
    out: Dict[str, Any] = {}
    for key, typ in keys.items():
        env_key = f"{prefix}_{key}"
        val = os.getenv(env_key)
        if val is None:
            continue
        if typ is int:
            try:
                out[key] = int(val)
            except ValueError:
                continue
        elif typ is float:
            try:
                out[key] = float(val)
            except ValueError:
                continue
        elif typ is bool:
            out[key] = _coerce_bool(val)
        else:
            out[key] = val
    return out


#
# Generic configuration loader used across the project
#


def _read_text(path: Path) -> str:
    with open(path, "rt", encoding="utf-8") as fh:
        return fh.read()


def load_mapping_from_file(config_path: str | Path) -> Dict[str, Any]:
    """Load a structured mapping from YAML, TOML or JSON file.

    - Supports .yaml/.yml (if PyYAML available), .toml (if tomllib available), and .json via stdlib
    - Returns a Python dict; raises ValueError if format unsupported or parsing fails
    """
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(f"Config file not found: {p}")

    suffix = p.suffix.lower()
    text = _read_text(p)

    if suffix in {".yaml", ".yml"}:
        if yaml is not None:  # Prefer PyYAML if available
            data = yaml.safe_load(text)
            if data is None:
                return {}
            if not isinstance(data, dict):
                raise ValueError("Top-level YAML must be a mapping")
            return dict(data)
        # Minimal fallback parser for our simple config files
        return _parse_simple_yaml_mapping(text)

    if suffix == ".toml":
        if tomllib is None:
            raise ValueError("TOML requested but tomllib unavailable. Use Python 3.11+ or use JSON/YAML.")
        return dict(tomllib.loads(text))

    if suffix == ".json":
        import json

        return dict(json.loads(text))

    raise ValueError(f"Unsupported config format: {suffix}")


def apply_env_overrides(config: Mapping[str, Any], *, prefix: str = "AK") -> Dict[str, Any]:
    """Apply simple environment overrides to a shallow config mapping.

    Supported keys via env:
    - {prefix}_THREADS -> int
    - {prefix}_WORK_DIR -> str
    - {prefix}_LOG_DIR -> str
    """
    out: Dict[str, Any] = dict(config)

    threads = os.getenv(f"{prefix}_THREADS")
    if threads is not None:
        try:
            out["threads"] = int(threads)
        except ValueError as e:
            # Log warning but continue with default
            import logging
            logging.getLogger(__name__).warning(f"Invalid threads value '{threads}', using default: {e}")

    work_dir = os.getenv(f"{prefix}_WORK_DIR")
    if work_dir:
        out["work_dir"] = work_dir

    log_dir = os.getenv(f"{prefix}_LOG_DIR")
    if log_dir:
        out["log_dir"] = log_dir

    return out


def _parse_simple_yaml_mapping(text: str) -> Dict[str, Any]:
    """Parse a minimal subset of YAML used in tests/configs when PyYAML is unavailable.

    Supported constructs:
      - key: value (ints and bare strings)
      - key: {} and key: []
      - inline dict: key: { a: 1, b: 2 }
      - simple nested mapping one or two levels deep with indentation
      - simple block lists under a key:
            key:
              - item1
              - item2
    This is not a general YAML parser.
    """
    result: Dict[str, Any] = {}
    lines = [ln.rstrip("\n") for ln in text.splitlines()]
    i = 0

    def indent_of(s: str) -> int:
        return len(s) - len(s.lstrip(" "))

    def parse_scalar(s: str) -> Any:
        s = s.strip()
        if s == "{}":
            return {}
        if s == "[]":
            return []
        try:
            return int(s)
        except ValueError:
            return s

    def parse_inline_dict(s: str) -> Dict[str, Any]:
        inner = s.strip()[1:-1].strip()
        if not inner:
            return {}
        items: Dict[str, Any] = {}
        parts = [p.strip() for p in inner.split(",") if p.strip()]
        for part in parts:
            if ":" not in part:
                continue
            k, v = part.split(":", 1)
            items[k.strip()] = parse_scalar(v)
        return items

    def collect_block_list(start_index: int, parent_indent: int) -> tuple[list[Any], int]:
        items: list[Any] = []
        j = start_index
        while j < len(lines):
            nxt = lines[j]
            if indent_of(nxt) <= parent_indent:
                break
            stripped = nxt.lstrip()
            if stripped.startswith("- "):
                items.append(parse_scalar(stripped[2:]))
                j += 1
                continue
            break
        return items, j

    def collect_child_mapping(start_index: int, parent_indent: int, depth: int = 1) -> tuple[Dict[str, Any], int]:
        children: Dict[str, Any] = {}
        j = start_index
        while j < len(lines):
            nxt = lines[j]
            if indent_of(nxt) <= parent_indent:
                break
            stripped = nxt.lstrip()
            if ":" not in stripped:
                j += 1
                continue
            ckey, crest = stripped.split(":", 1)
            ckey = ckey.strip()
            crest = crest.strip()
            if crest:
                if crest.startswith("{") and crest.endswith("}"):
                    children[ckey] = parse_inline_dict(crest)
                else:
                    children[ckey] = parse_scalar(crest)
                j += 1
                continue
            # Empty value: could be a list or a nested mapping one level deeper
            # Check for list items
            if j + 1 < len(lines) and lines[j + 1].lstrip().startswith("- "):
                items, j2 = collect_block_list(j + 1, indent_of(nxt))
                children[ckey] = items
                j = j2
                continue
            # Nested mapping one level deeper (only descend once more to keep simple)
            if depth < 2:
                grand_children, j2 = collect_child_mapping(j + 1, indent_of(nxt), depth=depth + 1)
                children[ckey] = grand_children
                j = j2
            else:
                children[ckey] = {}
                j += 1
        return children, j

    while i < len(lines):
        line = lines[i]
        i += 1
        if not line or line.lstrip().startswith("#"):
            continue
        if ":" not in line:
            continue
        key, rest = line.split(":", 1)
        key = key.strip()
        rest = rest.strip()
        if rest:
            if rest.startswith("{") and rest.endswith("}"):
                result[key] = parse_inline_dict(rest)
            else:
                result[key] = parse_scalar(rest)
            continue
        # No inline value: could be a list or mapping
        parent_indent = indent_of(line)
        # Block list under this key
        if i < len(lines) and lines[i].lstrip().startswith("- "):
            items, i = collect_block_list(i, parent_indent)
            result[key] = items
            continue
        # Child mapping (support up to two levels deep)
        children, i = collect_child_mapping(i, parent_indent, depth=1)
        result[key] = children

    return result
