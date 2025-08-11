from __future__ import annotations

import os
from typing import Any

from metainformant.core import config as core_config


def test_load_postgres_config_from_env_roundtrip() -> None:
    backup = {k: os.environ.get(k) for k in ["PG_HOST", "PG_DATABASE", "PG_USER", "PG_PASSWORD", "PG_PORT"]}
    os.environ["PG_HOST"] = "localhost"
    os.environ["PG_DATABASE"] = "db"
    os.environ["PG_USER"] = "user"
    os.environ["PG_PASSWORD"] = "pass"
    os.environ["PG_PORT"] = "5433"
    cfg = core_config.load_postgres_config_from_env()
    assert cfg and cfg.host == "localhost" and cfg.port == 5433
    for k, v in backup.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v


def test_load_typed_env() -> None:
    prefix = "APP"
    backup = {k: os.environ.get(k) for k in [f"{prefix}_STR", f"{prefix}_INT", f"{prefix}_BOOL_T", f"{prefix}_BOOL_F"]}
    os.environ[f"{prefix}_STR"] = "hello"
    os.environ[f"{prefix}_INT"] = "42"
    os.environ[f"{prefix}_BOOL_T"] = "true"
    os.environ[f"{prefix}_BOOL_F"] = "0"

    values = core_config.load_typed_env(prefix=prefix, keys={
        "STR": str,
        "INT": int,
        "BOOL_T": bool,
        "BOOL_F": bool,
    })

    assert values == {"STR": "hello", "INT": 42, "BOOL_T": True, "BOOL_F": False}
    for k, v in backup.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v


