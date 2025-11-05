"""Tests for core database utilities (optional PostgreSQL integration)."""

from __future__ import annotations

import os
import pytest

from metainformant.core.db import (
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
    """Test database connection functionality (requires PostgreSQL)."""

    @pytest.mark.skipif(
        os.environ.get("TEST_DATABASE_AVAILABLE") != "1",
        reason="PostgreSQL not available - real implementation requires database connection"
    )
    def test_get_db_client_with_config(self):
        """Test database client creation with environment configuration.
        
        Requires:
            - PG_HOST, PG_PORT, PG_DATABASE, PG_USER, PG_PASSWORD environment variables
            - OR: DB_NAME, DB_USER, DB_PASSWORD, PG_HOST environment variables
        """
        try:
            conn, cur = get_db_client()
            assert conn is not None
            assert cur is not None
            
            # Test basic query
            cur.execute("SELECT 1")
            result = cur.fetchone()
            assert result == (1,)
            
            # Cleanup
            cur.close()
            conn.close()
        except RuntimeError as e:
            if "Postgres configuration not found" in str(e):
                pytest.skip("PostgreSQL configuration not available in environment")
            raise

    def test_get_db_client_without_config(self):
        """Test database client creation fails gracefully without configuration."""
        # Save original environment
        original_env = {}
        for key in ["PG_HOST", "PG_PORT", "PG_DATABASE", "PG_USER", "PG_PASSWORD", 
                   "DB_NAME", "DB_USER", "DB_PASSWORD"]:
            original_env[key] = os.environ.get(key)
            if key in os.environ:
                del os.environ[key]
        
        try:
            with pytest.raises(RuntimeError, match="Postgres configuration not found"):
                get_db_client()
        finally:
            # Restore original environment
            for key, value in original_env.items():
                if value is not None:
                    os.environ[key] = value




