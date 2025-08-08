from __future__ import annotations

import os
from typing import Any

from metainformant.core import config as core_config


def test_load_postgres_config_from_env_roundtrip(monkeypatch: Any) -> None:
    monkeypatch.setenv("PG_HOST", "localhost")
    monkeypatch.setenv("PG_DATABASE", "db")
    monkeypatch.setenv("PG_USER", "user")
    monkeypatch.setenv("PG_PASSWORD", "pass")
    monkeypatch.setenv("PG_PORT", "5433")
    cfg = core_config.load_postgres_config_from_env()
    assert cfg and cfg.host == "localhost" and cfg.port == 5433


def test_load_typed_env(monkeypatch: Any) -> None:
    prefix = "APP"
    monkeypatch.setenv(f"{prefix}_STR", "hello")
    monkeypatch.setenv(f"{prefix}_INT", "42")
    monkeypatch.setenv(f"{prefix}_BOOL_T", "true")
    monkeypatch.setenv(f"{prefix}_BOOL_F", "0")

    values = core_config.load_typed_env(prefix=prefix, keys={
        "STR": str,
        "INT": int,
        "BOOL_T": bool,
        "BOOL_F": bool,
    })

    assert values == {"STR": "hello", "INT": 42, "BOOL_T": True, "BOOL_F": False}


