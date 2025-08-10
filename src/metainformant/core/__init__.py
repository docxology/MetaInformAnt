"""Shared cross-domain utilities for METAINFORMANT.

Modules:
- config: configuration helpers (env loaders, database config)
- db: database client helpers (optional)
- io: robust I/O utilities (json, jsonl, csv/tsv, gzip-aware)
- logging: consistently formatted loggers
- text: small text normalization helpers
- parallel: simple thread-based map preserving order
- hash: content and file hashing
"""

# Import core utilities eagerly
from . import config, io, logging, text, parallel, hash  # noqa: F401

# Import optional modules defensively to avoid hard deps during unrelated commands
try:  # noqa: SIM105
    from . import db  # type: ignore  # noqa: F401
except Exception:  # pragma: no cover - optional dependency may be missing
    db = None  # type: ignore


