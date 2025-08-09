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
        except ValueError:
            pass

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
      - simple nested mapping one level deep with indentation
    This is not a general YAML parser.
    """
    result: Dict[str, Any] = {}
    lines = [ln.rstrip("\n") for ln in text.splitlines()]
    i = 0

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
        # Collect indented child lines
        children: Dict[str, Any] = {}
        while i < len(lines):
            nxt = lines[i]
            if not nxt.startswith(" "):
                break
            stripped = nxt.lstrip()
            if ":" in stripped:
                ckey, crest = stripped.split(":", 1)
                ckey = ckey.strip()
                crest = crest.strip()
                if crest:
                    if crest.startswith("{") and crest.endswith("}"):
                        children[ckey] = parse_inline_dict(crest)
                    else:
                        children[ckey] = parse_scalar(crest)
                else:
                    children[ckey] = {}
            i += 1
        result[key] = children

    return result

