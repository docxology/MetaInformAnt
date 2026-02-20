"""Core Input/Output utilities for METAINFORMANT.

This module provides a unified interface for file operations, caching,
downloads, and disk management."""
from __future__ import annotations

from . import atomic, cache, checksums, disk, download, download_manager, download_robust, errors, io, paths

# Re-export public functions from io.py
from .io import (
    batch_download, download_csv, download_file, download_json, download_text,
    dump_json, dump_json_gz, dump_yaml,
    ensure_directory,
    load_json, load_json_gz, load_toml, load_yaml,
    open_text_auto,
    read_csv, read_delimited, read_jsonl, read_parquet, read_tsv,
    write_csv, write_delimited, write_jsonl, write_parquet, write_tsv,
)

# Re-export public functions from paths.py
from .paths import (
    change_extension, expand_and_resolve, find_files_by_extension,
    get_cache_dir, get_data_dir, get_directory_size, get_file_extension,
    get_file_size, get_logs_dir, get_project_root, get_temp_dir,
    is_safe_path, is_within, prepare_file_path, resolve_path,
    sanitize_filename, create_temp_file,
)

# Re-export download utilities
from .download import download_with_progress, DownloadResult

__all__ = [
    # Submodules
    'atomic', 'cache', 'checksums', 'disk', 'download', 'download_manager',
    'download_robust', 'errors', 'io', 'paths',
    # io.py functions
    'batch_download', 'download_csv', 'download_file', 'download_json', 'download_text',
    'dump_json', 'dump_json_gz', 'dump_yaml',
    'ensure_directory',
    'load_json', 'load_json_gz', 'load_toml', 'load_yaml',
    'open_text_auto',
    'read_csv', 'read_delimited', 'read_jsonl', 'read_parquet', 'read_tsv',
    'write_csv', 'write_delimited', 'write_jsonl', 'write_parquet', 'write_tsv',
    # paths.py functions
    'change_extension', 'expand_and_resolve', 'find_files_by_extension',
    'get_cache_dir', 'get_data_dir', 'get_directory_size', 'get_file_extension',
    'get_file_size', 'get_logs_dir', 'get_project_root', 'get_temp_dir',
    'is_safe_path', 'is_within', 'prepare_file_path', 'resolve_path',
    'sanitize_filename', 'create_temp_file',
    # download.py
    'download_with_progress', 'DownloadResult',
]
