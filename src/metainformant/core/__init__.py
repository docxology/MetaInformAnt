"""Shared cross-domain utilities for METAINFORMANT.

Modules:
- config: configuration helpers (env loaders, database config)
- db: database client helpers (optional)
- io: robust I/O utilities (json, jsonl, csv/tsv, gzip-aware)
- logging: consistently formatted loggers
- text: small text normalization helpers
- parallel: simple thread-based map preserving order
- hash: content and file hashing
- paths: path expansion, resolution, and containment checks
- cache: simple JSON cache with TTL
- errors: error handling and resilience utilities
- progress: progress tracking utilities
- validation: validation utilities
- workflow: workflow orchestration functions
"""

from __future__ import annotations

from pathlib import Path

# Import core utilities eagerly
from . import cache, config, disk, errors, hash, io, logging, parallel, paths, progress, text, validation, workflow  # noqa: F401

# Import optional modules defensively to avoid hard deps during unrelated commands
try:  # noqa: SIM105
    from . import db  # type: ignore  # noqa: F401
except Exception:  # pragma: no cover - optional dependency may be missing
    db = None  # type: ignore

# Type aliases for function signatures
PathType = str | Path | None

# Re-export workflow functions for backward compatibility
from .workflow import (  # noqa: F401
    create_sample_config,
    download_and_process_data,
    run_config_based_workflow,
    validate_config_file,
)
