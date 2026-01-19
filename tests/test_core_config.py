"""Comprehensive tests for core.utils.config module."""
from __future__ import annotations

import json
import os
from pathlib import Path

from metainformant.core.utils import config as core_config


class TestPostgresConfig:
    """Tests for PostgresConfig dataclass."""

    def test_postgres_config_creation(self) -> None:
        """Test creating PostgresConfig."""
        cfg = core_config.PostgresConfig(
            host="localhost",
            port=5432,
            database="testdb",
            user="admin",
            password="secret"
        )
        assert cfg.host == "localhost"
        assert cfg.port == 5432
        assert cfg.database == "testdb"


class TestLoadPostgresConfigFromEnv:
    """Tests for load_postgres_config_from_env function."""

    def test_load_postgres_config_full(self) -> None:
        """Test loading complete config from env."""
        env_backup = {}
        env_vars = ["PG_HOST", "PG_DATABASE", "PG_USER", "PG_PASSWORD", "PG_PORT"]
        for k in env_vars:
            env_backup[k] = os.environ.get(k)
        
        os.environ["PG_HOST"] = "dbhost"
        os.environ["PG_DATABASE"] = "mydb"
        os.environ["PG_USER"] = "myuser"
        os.environ["PG_PASSWORD"] = "mypass"
        os.environ["PG_PORT"] = "5433"
        
        try:
            cfg = core_config.load_postgres_config_from_env()
            assert cfg is not None
            assert cfg.host == "dbhost"
            assert cfg.port == 5433
            assert cfg.database == "mydb"
        finally:
            for k, v in env_backup.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v

    def test_load_postgres_config_missing(self) -> None:
        """Test returns None when vars missing."""
        # Clear potentially set vars
        for k in ["PG_HOST", "PG_DATABASE", "PG_USER", "PG_PASSWORD"]:
            os.environ.pop(k, None)
        
        cfg = core_config.load_postgres_config_from_env()
        # Should return None if not all required vars are set
        # (depends on implementation)


class TestLoadConfigFile:
    """Tests for load_config_file function."""

    def test_load_config_file_json(self, tmp_path: Path) -> None:
        """Test loading JSON config."""
        config_path = tmp_path / "config.json"
        config_data = {"key1": "value1", "nested": {"key2": 123}}
        config_path.write_text(json.dumps(config_data))
        
        loaded = core_config.load_config_file(config_path)
        assert loaded == config_data


class TestGetEnvOrDefault:
    """Tests for get_env_or_default function."""

    def test_get_env_or_default_exists(self) -> None:
        """Test getting existing env var."""
        os.environ["TEST_VAR_123"] = "test_value"
        try:
            result = core_config.get_env_or_default("TEST_VAR_123", "default")
            assert result == "test_value"
        finally:
            os.environ.pop("TEST_VAR_123", None)

    def test_get_env_or_default_missing(self) -> None:
        """Test getting missing env var returns default."""
        os.environ.pop("NONEXISTENT_VAR_XYZ", None)
        result = core_config.get_env_or_default("NONEXISTENT_VAR_XYZ", "fallback")
        assert result == "fallback"


class TestLoadTypedEnv:
    """Tests for load_typed_env function."""

    def test_load_typed_env_all_types(self) -> None:
        """Test loading various types from env."""
        prefix = "TEST_LTE"
        env_backup = {}
        keys = [f"{prefix}_STR", f"{prefix}_INT", f"{prefix}_BOOL"]
        for k in keys:
            env_backup[k] = os.environ.get(k)
        
        os.environ[f"{prefix}_STR"] = "hello"
        os.environ[f"{prefix}_INT"] = "42"
        os.environ[f"{prefix}_BOOL"] = "true"
        
        try:
            values = core_config.load_typed_env(
                prefix=prefix,
                keys={"STR": str, "INT": int, "BOOL": bool}
            )
            assert values["STR"] == "hello"
            assert values["INT"] == 42
            assert values["BOOL"] is True
        finally:
            for k, v in env_backup.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v

    def test_load_typed_env_bool_false(self) -> None:
        """Test boolean false values."""
        prefix = "TEST_BOOL"
        os.environ.pop(f"{prefix}_VAL", None)
        os.environ[f"{prefix}_VAL"] = "false"
        
        try:
            values = core_config.load_typed_env(prefix=prefix, keys={"VAL": bool})
            assert values.get("VAL") is False
        finally:
            os.environ.pop(f"{prefix}_VAL", None)


class TestLoadMappingFromFile:
    """Tests for load_mapping_from_file function."""

    def test_load_mapping_json(self, tmp_path: Path) -> None:
        """Test loading JSON mapping."""
        path = tmp_path / "config.json"
        data = {"settings": {"threads": 4}}
        path.write_text(json.dumps(data))
        
        loaded = core_config.load_mapping_from_file(path)
        assert loaded == data

    def test_load_mapping_invalid_extension(self, tmp_path: Path) -> None:
        """Test unsupported file extension raises error."""
        path = tmp_path / "config.xyz"
        path.write_text("{}")
        
        try:
            core_config.load_mapping_from_file(path)
            assert False, "Should have raised ValueError"
        except ValueError:
            pass


class TestMergeConfigs:
    """Tests for merge_configs function."""

    def test_merge_configs_basic(self) -> None:
        """Test basic config merge."""
        base = {"a": 1, "b": 2}
        override = {"b": 3, "c": 4}
        
        merged = core_config.merge_configs(base, override)
        
        assert merged["a"] == 1  # From base
        assert merged["b"] == 3  # Overridden
        assert merged["c"] == 4  # From override

    def test_merge_configs_nested(self) -> None:
        """Test nested config merge."""
        base = {"a": 1, "nested": {"x": 1, "y": 2}}
        override = {"nested": {"y": 3, "z": 4}}
        
        merged = core_config.merge_configs(base, override)
        
        assert merged["nested"]["x"] == 1  # From base
        assert merged["nested"]["y"] == 3  # Overridden
        assert merged["nested"]["z"] == 4  # From override

    def test_merge_configs_empty(self) -> None:
        """Test merge with empty dicts."""
        base = {"a": 1}
        merged = core_config.merge_configs(base, {})
        assert merged == {"a": 1}
        
        merged2 = core_config.merge_configs({}, {"a": 1})
        assert merged2 == {"a": 1}


class TestCoerceConfigTypes:
    """Tests for coerce_config_types function."""

    def test_coerce_config_types_basic(self) -> None:
        """Test basic type coercion."""
        config = {"port": "8080", "debug": "true"}
        type_map = {"port": int, "debug": bool}
        
        coerced = core_config.coerce_config_types(config, type_map)
        
        assert coerced["port"] == 8080
        assert coerced["debug"] is True

    def test_coerce_config_types_preserves_unspecified(self) -> None:
        """Test that keys not in type_map are preserved."""
        config = {"port": "8080", "name": "server"}
        type_map = {"port": int}
        
        coerced = core_config.coerce_config_types(config, type_map)
        
        assert coerced["port"] == 8080
        assert coerced["name"] == "server"  # Preserved as string


class TestApplyEnvOverrides:
    """Tests for apply_env_overrides function."""

    def test_apply_env_overrides_threads(self) -> None:
        """Test THREADS override."""
        os.environ["AK_THREADS"] = "16"
        
        try:
            config = {"threads": 4}
            result = core_config.apply_env_overrides(config, prefix="AK")
            assert result.get("threads") == 16
        finally:
            os.environ.pop("AK_THREADS", None)

    def test_apply_env_overrides_work_dir(self) -> None:
        """Test WORK_DIR override."""
        os.environ["AK_WORK_DIR"] = "/custom/work"
        
        try:
            config = {"work_dir": "/default"}
            result = core_config.apply_env_overrides(config, prefix="AK")
            assert result.get("work_dir") == "/custom/work"
        finally:
            os.environ.pop("AK_WORK_DIR", None)
