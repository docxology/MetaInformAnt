"""Core Input/Output utilities for METAINFORMANT.

This module provides a unified interface for file operations, caching,
downloads, and disk management.
"""

from .errors import (
    IOError,
    FileNotFoundError,
    CacheError,
    DownloadError,
)

from .io import (
    ensure_directory,
    open_text_auto,
    load_json,
    dump_json,
    dump_json_gz,
    load_json_gz,
    load_yaml,
    load_toml,
    dump_yaml,
    read_parquet,
    write_parquet,
    read_jsonl,
    write_jsonl,
    read_delimited,
    write_delimited,
    read_csv,
    write_csv,
    read_tsv,
    write_tsv,
    download_file,
    download_json,
    download_text,
    download_csv,
    batch_download,
)

from .cache import (
    JsonCache,
    CacheEntry,
    get_cache_info,
    cache_json,
    load_cached_json,
    clear_cache_dir,
)

from .disk import (
    get_disk_usage,
    get_free_space,
    get_recommended_temp_dir,
    ensure_disk_space,
    check_disk_space,
    cleanup_temp_files,
    cleanup_old_files,
    get_directory_size,
    get_largest_files,
    monitor_disk_space,
    safe_remove_directory,
    get_disk_space_info,
    detect_drive_size_category,
    get_recommended_batch_size,
)

from .paths import (
    get_project_root,
    get_data_dir,
    get_cache_dir,
    get_logs_dir,
    get_temp_dir,
    resolve_path,
    is_safe_path,
)

__all__ = [
    # Errors
    "IOError",
    "FileNotFoundError",
    "CacheError",
    "DownloadError",
    # IO
    "ensure_directory",
    "open_text_auto",
    "load_json",
    "dump_json",
    "dump_json_gz",
    "load_json_gz",
    "load_yaml",
    "load_toml",
    "dump_yaml",
    "read_parquet",
    "write_parquet",
    "read_jsonl",
    "write_jsonl",
    "read_delimited",
    "write_delimited",
    "read_csv",
    "write_csv",
    "read_tsv",
    "write_tsv",
    "download_file",
    "download_json",
    "download_text",
    "download_csv",
    "batch_download",
    # Cache
    "JsonCache",
    "CacheEntry",
    "get_cache_info",
    "cache_json",
    "load_cached_json",
    "clear_cache_dir",
    # Disk
    "get_disk_usage",
    "get_free_space",
    "get_recommended_temp_dir",
    "ensure_disk_space",
    "check_disk_space",
    "cleanup_temp_files",
    "cleanup_old_files",
    "get_directory_size",
    "get_largest_files",
    "monitor_disk_space",
    "safe_remove_directory",
    "get_disk_space_info",
    "detect_drive_size_category",
    "get_recommended_batch_size",
    # Paths
    "get_project_root",
    "get_data_dir",
    "get_cache_dir",
    "get_logs_dir",
    "get_temp_dir",
    "resolve_path",
    "is_safe_path",
]
