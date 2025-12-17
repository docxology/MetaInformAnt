from __future__ import annotations

from typing import Tuple

try:
    import psycopg2
    HAS_PSYCOPG2 = True
except ImportError:
    psycopg2 = None  # type: ignore
    HAS_PSYCOPG2 = False

from .config import load_postgres_config_from_env
from .logging import get_logger

logger = get_logger(__name__)


def get_db_client() -> Tuple["psycopg2.extensions.connection", "psycopg2.extensions.cursor"]:
    """Get PostgreSQL database client connection and cursor.

    Loads configuration from environment variables and establishes connection.

    Returns:
        Tuple of (connection, cursor) objects

    Raises:
        RuntimeError: If PostgreSQL configuration not found in environment
        ImportError: If psycopg2 is not available
    """
    if not HAS_PSYCOPG2:
        raise ImportError("psycopg2 is required for PostgreSQL database operations. Install with: uv add psycopg2-binary")

    config = load_postgres_config_from_env()
    if config is None:
        raise RuntimeError("Postgres configuration not found in environment variables")

    logger.debug("Initializing connection to PostgreSQL database")
    conn = psycopg2.connect(
        host=config.host,
        port=config.port,
        database=config.database,
        user=config.user,
        password=config.password,
    )
    cur = conn.cursor()
    cur.execute("SELECT 1")
    return conn, cur


def build_postgres_url(host: str, port: int, database: str, user: str, password: str) -> str:
    """Build PostgreSQL connection URL.

    Args:
        host: Database host
        port: Database port
        database: Database name
        user: Username
        password: Password

    Returns:
        PostgreSQL connection URL
    """
    return f"postgresql://{user}:{password}@{host}:{port}/{database}"


def sanitize_connection_params(params: dict) -> dict:
    """Sanitize database connection parameters for security.

    Args:
        params: Raw connection parameters

    Returns:
        Sanitized connection parameters
    """
    import re

    sanitized = {}
    for key, value in params.items():
        if isinstance(value, str):
            # Remove SQL injection attempts and dangerous characters
            cleaned = re.sub(r'[;\'"\\]', "", value)
            cleaned = re.sub(r"\b(DROP|DELETE|INSERT|UPDATE|CREATE|ALTER)\b", "", cleaned, flags=re.IGNORECASE)
            cleaned = cleaned.strip()
            sanitized[key] = cleaned
        else:
            sanitized[key] = value

    return sanitized
