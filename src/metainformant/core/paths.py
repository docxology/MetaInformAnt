"""Deprecated compatibility facade for path helpers.

Use :mod:`metainformant.core.io.paths` for new code.  This module keeps the
historic ``metainformant.core.paths`` import path working while documentation,
scripts, and downstream callers migrate to the canonical location.
"""

from __future__ import annotations

from metainformant.core.io.paths import *  # noqa: F403

__all__ = [
    "expand_and_resolve",
    "is_within",
    "ensure_directory",
    "prepare_file_path",
    "is_safe_path",
    "get_file_extension",
    "change_extension",
    "find_files_by_extension",
    "get_file_size",
    "get_directory_size",
    "sanitize_filename",
    "create_temp_file",
    "discover_output_patterns",
    "find_output_locations",
    "get_module_output_base",
    "list_output_structure",
    "get_project_root",
    "get_data_dir",
    "get_cache_dir",
    "get_logs_dir",
    "get_temp_dir",
    "resolve_path",
]
