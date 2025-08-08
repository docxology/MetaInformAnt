from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Callable, Dict, Mapping


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

