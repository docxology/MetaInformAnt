"""Deprecated compatibility facade for configuration helpers.

Use :mod:`metainformant.core.utils.config` for new code.  This module keeps the
historic ``metainformant.core.config`` import path working while the repository
finishes migrating examples and downstream callers.
"""

from __future__ import annotations

from metainformant.core.utils.config import *  # noqa: F403

__all__ = [
    "PostgresConfig",
    "load_postgres_config_from_env",
    "load_config_file",
    "get_env_or_default",
    "load_typed_env",
    "load_mapping_from_file",
    "merge_configs",
    "coerce_config_types",
    "apply_env_overrides",
    "discover_config_files",
    "get_config_schema",
    "find_configs_for_module",
    "list_config_templates",
]
