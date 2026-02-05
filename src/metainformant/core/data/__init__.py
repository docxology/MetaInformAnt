"""Core data utilities for METAINFORMANT.

This module provides data validation, database integration, and type checking utilities.
"""

from __future__ import annotations

# Database utilities (optional - requires psycopg2)
try:
    from .db import (
        PostgresConnection,
        build_postgres_url,
        get_connection,
        get_db_client,
        sanitize_connection_params,
    )

    _HAS_DB = True
except ImportError:
    _HAS_DB = False
    PostgresConnection = None
    build_postgres_url = None
    get_db_client = None
    sanitize_connection_params = None
    get_connection = None

# Validation utilities
from .validation import (
    validate_json_schema,
    validate_not_empty,
    validate_not_none,
    validate_path_exists,
    validate_path_is_dir,
    validate_path_is_file,
    validate_path_within,
    validate_range,
    validate_schema,
    validate_type,
    validator,
)

__all__ = [
    # Database (optional)
    "PostgresConnection",
    "build_postgres_url",
    "get_db_client",
    "sanitize_connection_params",
    "get_connection",
    # Validation
    "validate_type",
    "validate_range",
    "validate_path_exists",
    "validate_path_is_file",
    "validate_path_is_dir",
    "validate_path_within",
    "validate_not_none",
    "validate_not_empty",
    "validate_schema",
    "validate_json_schema",
    "validator",
]
