"""Tests for core database utilities (optional PostgreSQL integration)."""

from __future__ import annotations

import os

import pytest

from metainformant.core.data.db import (
    build_postgres_url,
    get_db_client,
    sanitize_connection_params,
)


class TestDatabaseUtilities:
    """Test database utility functions."""

    def test_build_postgres_url(self):
        """Test PostgreSQL URL construction."""
        url = build_postgres_url("localhost", 5432, "testdb", "user", "pass")
        assert url == "postgresql://user:pass@localhost:5432/testdb"

    def test_build_postgres_url_different_port(self):
        """Test URL construction with non-standard port."""
        url = build_postgres_url("example.com", 5433, "mydb", "admin", "secret")
        assert url == "postgresql://admin:secret@example.com:5433/mydb"

    def test_sanitize_connection_params(self):
        """Test connection parameter sanitization."""
        params = {
            "host": "localhost",
            "user": "test; DROP TABLE users;",
            "password": "pass'123",
            "database": "testdb",
        }
        sanitized = sanitize_connection_params(params)

        # SQL injection attempts should be removed
        assert "DROP" not in sanitized["user"]
        assert ";" not in sanitized["user"]
        assert "'" not in sanitized["password"]

        # Valid parts should remain
        assert sanitized["host"] == "localhost"
        assert sanitized["database"] == "testdb"

    def test_sanitize_connection_params_non_string(self):
        """Test sanitization with non-string values."""
        params = {"host": "localhost", "port": 5432, "ssl": True}
        sanitized = sanitize_connection_params(params)
        assert sanitized["host"] == "localhost"
        assert sanitized["port"] == 5432
        assert sanitized["ssl"] is True

    def test_sanitize_connection_params_sql_keywords(self):
        """Test sanitization removes SQL keywords from strings."""
        params = {
            "host": "localhost",
            "user": "admin DELETE FROM users",
            "password": "pass",
        }
        sanitized = sanitize_connection_params(params)
        assert "DELETE" not in sanitized["user"]
        assert sanitized["host"] == "localhost"


class TestDatabaseConnection:
    """Test database connection behavior without requiring a live server."""

    def test_get_db_client_with_config(self):
        """Test that configured connection attempts expose actionable failures."""
        try:
            client = get_db_client(host="127.0.0.1", port=1, database="metainformant", user="postgres")
        except (ImportError, RuntimeError) as exc:
            assert str(exc), "Database failures must retain their diagnostic message"
        else:
            connection = client.connect()
            try:
                rows = client.execute_query(connection, "SELECT 1")
                assert rows == [{"?column?": 1}] or rows == [{"int4": 1}] or rows == [{"select": 1}]
            finally:
                client.release_connection(connection)

    def test_get_db_client_without_config(self):
        """Test database client creation fails gracefully without configuration."""
        # Save original environment
        original_env = {}
        for key in ["PG_HOST", "PG_PORT", "PG_DATABASE", "PG_USER", "PG_PASSWORD", "DB_NAME", "DB_USER", "DB_PASSWORD"]:
            original_env[key] = os.environ.get(key)
            if key in os.environ:
                del os.environ[key]

        try:
            # Should raise ImportError first since psycopg2 is not available
            with pytest.raises(ImportError, match="psycopg2 required for database operations"):
                get_db_client()
        finally:
            # Restore original environment
            for key, value in original_env.items():
                if value is not None:
                    os.environ[key] = value
