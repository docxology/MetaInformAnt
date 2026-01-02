"""Database integration helpers for METAINFORMANT.

This module provides PostgreSQL database connection pooling, query execution,
and transaction management for bioinformatics data storage and retrieval.
"""

from __future__ import annotations

import contextlib
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union

from metainformant.core.logging import get_logger

logger = get_logger(__name__)

# Optional imports for PostgreSQL support
try:
    import psycopg2
    import psycopg2.extras
    import psycopg2.pool
    HAS_PSYCOPG2 = True
except ImportError:
    psycopg2 = None
    HAS_PSYCOPG2 = False
    logger.warning("psycopg2 not available, PostgreSQL functionality disabled")


class PostgresConnection:
    """PostgreSQL database connection wrapper with connection pooling.

    This class provides a high-level interface for PostgreSQL database operations
    with automatic connection management and error handling.

    Example:
        # Create connection
        conn = PostgresConnection(
            host="localhost",
            port=5432,
            database="metainformant",
            user="user",
            password="password"
        )

        # Use connection
        with conn.connect() as db_conn:
            results = conn.execute_query(
                db_conn,
                "SELECT * FROM genes WHERE species = %s",
                ("human",)
            )
    """

    def __init__(
        self,
        host: str = "localhost",
        port: int = 5432,
        database: str = "metainformant",
        user: str = "postgres",
        password: str = "",
        min_connections: int = 1,
        max_connections: int = 10,
    ):
        """Initialize PostgreSQL connection pool.

        Args:
            host: Database host
            port: Database port
            database: Database name
            user: Database user
            password: Database password
            min_connections: Minimum pool size
            max_connections: Maximum pool size
        """
        if not HAS_PSYCOPG2:
            raise ImportError("psycopg2 required for PostgreSQL functionality")

        self.host = host
        self.port = port
        self.database = database
        self.user = user
        self.password = password

        # Connection pool
        self._pool = psycopg2.pool.ThreadedConnectionPool(
            minconn=min_connections,
            maxconn=max_connections,
            host=host,
            port=port,
            database=database,
            user=user,
            password=password,
        )

        logger.info(f"Initialized PostgreSQL connection pool to {host}:{port}/{database}")

    def connect(self) -> psycopg2.extensions.connection:
        """Get connection from pool.

        Returns:
            Database connection

        Raises:
            RuntimeError: If connection pool not available
        """
        try:
            return self._pool.getconn()
        except Exception as e:
            logger.error(f"Failed to get database connection: {e}")
            raise RuntimeError(f"Database connection failed: {e}") from e

    def release_connection(self, conn: psycopg2.extensions.connection) -> None:
        """Release connection back to pool.

        Args:
            conn: Connection to release
        """
        try:
            self._pool.putconn(conn)
        except Exception as e:
            logger.warning(f"Failed to release connection: {e}")

    def execute_query(
        self,
        conn: psycopg2.extensions.connection,
        query: str,
        params: Tuple = (),
        fetch: bool = True,
    ) -> List[Dict[str, Any]]:
        """Execute SQL query.

        Args:
            conn: Database connection
            query: SQL query string
            params: Query parameters
            fetch: Whether to fetch results

        Returns:
            List of result dictionaries if fetch=True, empty list otherwise

        Raises:
            RuntimeError: If query execution fails
        """
        try:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:
                cursor.execute(query, params)

                if fetch and cursor.description:
                    return [dict(row) for row in cursor.fetchall()]
                else:
                    return []

        except Exception as e:
            logger.error(f"Query execution failed: {query} with params {params}: {e}")
            conn.rollback()
            raise RuntimeError(f"Database query failed: {e}") from e

    def execute_transaction(
        self,
        conn: psycopg2.extensions.connection,
        queries: List[Tuple[str, Tuple]],
    ) -> List[List[Dict[str, Any]]]:
        """Execute multiple queries in a transaction.

        Args:
            conn: Database connection
            queries: List of (query, params) tuples

        Returns:
            List of result lists for each query

        Raises:
            RuntimeError: If transaction fails
        """
        results = []
        try:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:
                for query, params in queries:
                    cursor.execute(query, params)
                    if cursor.description:
                        results.append([dict(row) for row in cursor.fetchall()])
                    else:
                        results.append([])

            conn.commit()
            return results

        except Exception as e:
            conn.rollback()
            logger.error(f"Transaction failed: {e}")
            raise RuntimeError(f"Database transaction failed: {e}") from e

    def bulk_insert(
        self,
        conn: psycopg2.extensions.connection,
        table: str,
        columns: List[str],
        data: List[Tuple],
        batch_size: int = 1000,
    ) -> int:
        """Bulk insert data into table.

        Args:
            conn: Database connection
            table: Table name
            columns: Column names
            data: List of data tuples
            batch_size: Batch size for insertion

        Returns:
            Number of rows inserted

        Raises:
            RuntimeError: If bulk insert fails
        """
        if not data:
            return 0

        try:
            with conn.cursor() as cursor:
                # Create placeholders for bulk insert
                placeholders = ", ".join(["%s"] * len(columns))
                query = f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({placeholders})"

                # Insert in batches
                total_inserted = 0
                for i in range(0, len(data), batch_size):
                    batch = data[i : i + batch_size]
                    cursor.executemany(query, batch)
                    total_inserted += len(batch)

                conn.commit()
                return total_inserted

        except Exception as e:
            conn.rollback()
            logger.error(f"Bulk insert failed for table {table}: {e}")
            raise RuntimeError(f"Database bulk insert failed: {e}") from e

    def create_table(
        self,
        conn: psycopg2.extensions.connection,
        table_name: str,
        schema: Dict[str, str],
        indexes: Optional[List[str]] = None,
        if_not_exists: bool = True,
    ) -> None:
        """Create database table.

        Args:
            conn: Database connection
            table_name: Name of table to create
            schema: Dictionary mapping column names to SQL types
            indexes: List of column names to index
            if_not_exists: Whether to use IF NOT EXISTS

        Raises:
            RuntimeError: If table creation fails
        """
        try:
            exists_clause = "IF NOT EXISTS" if if_not_exists else ""

            columns_def = ", ".join(f"{col} {col_type}" for col, col_type in schema.items())
            query = f"CREATE TABLE {exists_clause} {table_name} ({columns_def})"

            self.execute_query(conn, query, fetch=False)

            # Create indexes
            if indexes:
                for column in indexes:
                    if column in schema:
                        idx_name = f"idx_{table_name}_{column}"
                        idx_query = f"CREATE INDEX IF NOT EXISTS {idx_name} ON {table_name} ({column})"
                        self.execute_query(conn, idx_query, fetch=False)

            logger.info(f"Created table {table_name} with columns: {list(schema.keys())}")

        except Exception as e:
            logger.error(f"Failed to create table {table_name}: {e}")
            raise RuntimeError(f"Table creation failed: {e}") from e

    def table_exists(
        self,
        conn: psycopg2.extensions.connection,
        table_name: str,
    ) -> bool:
        """Check if table exists.

        Args:
            conn: Database connection
            table_name: Table name to check

        Returns:
            True if table exists
        """
        try:
            result = self.execute_query(
                conn,
                "SELECT EXISTS (SELECT 1 FROM information_schema.tables WHERE table_name = %s)",
                (table_name,),
            )
            return result[0]["exists"] if result else False

        except Exception as e:
            logger.warning(f"Failed to check table existence: {e}")
            return False

    def close(self) -> None:
        """Close connection pool."""
        if hasattr(self, '_pool'):
            self._pool.closeall()
            logger.info("Closed PostgreSQL connection pool")


def build_postgres_url(
    host: str,
    port: int,
    database: str,
    user: str,
    password: str,
) -> str:
    """Build PostgreSQL connection URL.

    Args:
        host: Database host
        port: Database port
        database: Database name
        user: Database user
        password: Database password

    Returns:
        PostgreSQL connection URL string

    Example:
        >>> build_postgres_url("localhost", 5432, "mydb", "user", "pass")
        'postgresql://user:pass@localhost:5432/mydb'
    """
    return f"postgresql://{user}:{password}@{host}:{port}/{database}"


def get_db_client(
    host: str = "localhost",
    port: int = 5432,
    database: str = "metainformant",
    user: str = "postgres",
    password: str = "",
    **kwargs: Any,
) -> PostgresConnection:
    """Get database client instance.

    Args:
        host: Database host
        port: Database port
        database: Database name
        user: Database user
        password: Database password
        **kwargs: Additional connection parameters

    Returns:
        PostgresConnection instance

    Raises:
        ImportError: If psycopg2 is not available
    """
    if not HAS_PSYCOPG2:
        raise ImportError("psycopg2 required for database operations")

    return PostgresConnection(
        host=host,
        port=port,
        database=database,
        user=user,
        password=password,
        **kwargs
    )


def sanitize_connection_params(params: Dict[str, Any]) -> Dict[str, Any]:
    """Sanitize database connection parameters to prevent SQL injection.

    Args:
        params: Dictionary of connection parameters

    Returns:
        Sanitized parameter dictionary

    Example:
        >>> params = {"user": "test; DROP TABLE", "host": "localhost"}
        >>> sanitize_connection_params(params)
        {'user': 'test', 'host': 'localhost'}
    """
    sanitized = {}

    for key, value in params.items():
        if isinstance(value, str):
            # Remove potentially dangerous SQL characters
            clean_value = value.replace("'", "").replace(";", "").replace("--", "")
            # Remove SQL keywords that could be used for injection
            dangerous_patterns = [
                "DROP", "DELETE", "UPDATE", "INSERT", "ALTER",
                "CREATE", "TRUNCATE", "EXEC", "EXECUTE"
            ]
            for pattern in dangerous_patterns:
                clean_value = clean_value.replace(pattern, "")
            sanitized[key] = clean_value
        else:
            # Keep non-string values as-is
            sanitized[key] = value

    return sanitized


@contextlib.contextmanager
def get_connection(
    host: str = "localhost",
    port: int = 5432,
    database: str = "metainformant",
    user: str = "postgres",
    password: str = "",
    **kwargs: Any,
) -> Iterator[PostgresConnection]:
    """Context manager for PostgreSQL connections.

    Args:
        host: Database host
        port: Database port
        database: Database name
        user: Database user
        password: Database password
        **kwargs: Additional connection parameters

    Yields:
        PostgresConnection instance

    Example:
        with get_connection(database="mydb") as conn:
            with conn.connect() as db_conn:
                results = conn.execute_query(db_conn, "SELECT * FROM table")
    """
    conn = None
    try:
        conn = PostgresConnection(
            host=host,
            port=port,
            database=database,
            user=user,
            password=password,
            **kwargs
        )
        yield conn
    finally:
        if conn:
            conn.close()

