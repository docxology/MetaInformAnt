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
- discovery: symbolic mapping and context discovery utilities
- symbols: symbol indexing and cross-referencing utilities
"""

from __future__ import annotations

from pathlib import Path

# Import core utilities eagerly
from . import (  # noqa: F401
    cache,
    config,
    discovery,
    disk,
    errors,
    filesystem,
    hash,
    io,
    logging,
    parallel,
    paths,
    progress,
    symbols,
    text,
    validation,
    workflow,
)

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

# Re-export discovery functions
from .discovery import (  # noqa: F401
    build_call_graph,
    discover_configs,
    discover_functions,
    discover_output_patterns,
    discover_workflows,
    find_symbol_usage,
    get_module_dependencies,
)

# Re-export symbol functions
from .symbols import (  # noqa: F401
    find_symbol,
    find_symbol_references,
    fuzzy_find_symbol,
    get_symbol_metadata,
    get_symbol_signature,
    index_classes,
    index_functions,
)

# Re-export config discovery functions
from .config import (  # noqa: F401
    discover_config_files,
    find_configs_for_module,
    get_config_schema,
    list_config_templates,
)

# Re-export path discovery functions
from .paths import (  # noqa: F401
    discover_output_patterns as discover_path_output_patterns,
    find_output_locations,
    get_module_output_base,
    list_output_structure,
)
