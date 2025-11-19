"""Test that UV environment variables are automatically configured.

This test verifies that conftest.py automatically sets UV_CACHE_DIR
and UV_PROJECT_ENVIRONMENT based on filesystem type.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest


def test_uv_env_vars_configured():
    """Test that UV environment variables are set by conftest."""
    # These should be set by setup_test_environment fixture
    assert "UV_CACHE_DIR" in os.environ, "UV_CACHE_DIR should be set by conftest"
    assert "UV_PROJECT_ENVIRONMENT" in os.environ, "UV_PROJECT_ENVIRONMENT should be set by conftest"
    
    # Verify paths exist or can be created
    cache_dir = Path(os.environ["UV_CACHE_DIR"])
    venv_dir = Path(os.environ["UV_PROJECT_ENVIRONMENT"])
    
    # Cache directory should exist (created by get_uv_cache_dir)
    assert cache_dir.exists() or cache_dir.parent.exists(), f"UV cache directory should exist: {cache_dir}"
    
    # Verify paths are absolute
    assert cache_dir.is_absolute(), f"UV_CACHE_DIR should be absolute: {cache_dir}"
    assert venv_dir.is_absolute(), f"UV_PROJECT_ENVIRONMENT should be absolute: {venv_dir}"


def test_uv_env_respects_user_override(monkeypatch, tmp_path):
    """Test that user-set environment variables are respected."""
    # Set custom values using tmp_path (writable location)
    custom_cache = str(tmp_path / "custom-uv-cache")
    custom_venv = str(tmp_path / "custom-venv")
    
    monkeypatch.setenv("UV_CACHE_DIR", custom_cache)
    monkeypatch.setenv("UV_PROJECT_ENVIRONMENT", custom_venv)
    
    # Verify the filesystem utilities respect env vars
    from metainformant.core.filesystem import get_uv_cache_dir, get_venv_location
    
    repo_root = Path(__file__).parent.parent
    cache_dir = get_uv_cache_dir(repo_root)
    venv_dir = get_venv_location(repo_root)
    
    assert str(cache_dir) == custom_cache, f"Should respect user-set UV_CACHE_DIR, got {cache_dir}"
    assert str(venv_dir) == custom_venv, f"Should respect user-set UV_PROJECT_ENVIRONMENT, got {venv_dir}"

