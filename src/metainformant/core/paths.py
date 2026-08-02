"""Deprecated compatibility facade for path helpers.

Use :mod:`metainformant.core.io.paths` for new code.  This module keeps the
historic ``metainformant.core.paths`` import path working while documentation,
scripts, and downstream callers migrate to the canonical location.
"""

from __future__ import annotations

from metainformant.core.io.paths import (
    change_extension,
    create_temp_file,
    discover_output_patterns,
    ensure_directory,
    expand_and_resolve,
    find_files_by_extension,
    find_output_locations,
    get_cache_dir,
    get_data_dir,
    get_directory_size,
    get_file_extension,
    get_file_size,
    get_logs_dir,
    get_module_output_base,
    get_project_root,
    get_temp_dir,
    is_safe_path,
    is_within,
    list_output_structure,
    prepare_file_path,
    resolve_path,
    sanitize_filename,
)

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
