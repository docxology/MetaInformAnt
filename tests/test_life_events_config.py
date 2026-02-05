"""Tests for life_events configuration management."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from metainformant.core.io.io import dump_json
from metainformant.life_events import LifeEventsWorkflowConfig, load_life_events_config


def test_load_life_events_config_basic(tmp_path: Path) -> None:
    """Test loading minimal config from YAML."""
    config_file = tmp_path / "test_config.yaml"
    config_content = """
work_dir: output/life_events/test
threads: 4
"""
    config_file.write_text(config_content)

    config = load_life_events_config(config_file)
    assert isinstance(config, LifeEventsWorkflowConfig)
    assert config.work_dir == Path("output/life_events/test").expanduser().resolve()
    assert config.threads == 4
    assert config.log_dir is None
    assert config.embedding == {}
    assert config.model == {}
    assert config.workflow == {}
    assert config.output == {}


def test_load_life_events_config_all_sections(tmp_path: Path) -> None:
    """Test loading complete config with all sections."""
    config_file = tmp_path / "full_config.yaml"
    config_content = """
work_dir: output/life_events/full_test
threads: 8
log_dir: output/life_events/logs

embedding:
  embedding_dim: 100
  window_size: 5
  epochs: 10
  method: skipgram
  learning_rate: 0.01

model:
  model_type: embedding
  task_type: classification
  random_state: 42

workflow:
  save_model: true
  save_embeddings: true

output:
  format: json
  include_probabilities: true
"""
    config_file.write_text(config_content)

    config = load_life_events_config(config_file)
    assert config.threads == 8
    assert config.log_dir is not None
    assert "embedding_dim" in config.embedding
    assert config.embedding["embedding_dim"] == 100
    assert config.embedding["window_size"] == 5
    assert config.model["model_type"] == "embedding"
    assert config.model["task_type"] == "classification"
    assert config.workflow["save_model"] is True
    assert config.output["format"] == "json"


def test_load_life_events_config_json(tmp_path: Path) -> None:
    """Test loading config from JSON format."""
    config_file = tmp_path / "test_config.json"
    config_data = {
        "work_dir": "output/life_events/test",
        "threads": 6,
        "embedding": {"embedding_dim": 100, "window_size": 5},
        "model": {"model_type": "embedding", "random_state": 42},
    }
    dump_json(config_data, config_file)

    config = load_life_events_config(config_file)
    assert config.threads == 6
    assert config.embedding["embedding_dim"] == 100
    assert config.model["model_type"] == "embedding"
    assert config.model["random_state"] == 42


def test_load_life_events_config_toml(tmp_path: Path) -> None:
    """Test loading config from TOML format (if available)."""
    try:
        import tomllib
    except ImportError:
        pytest.skip("TOML support not available (requires Python 3.11+)")

    config_file = tmp_path / "test_config.toml"
    config_content = """
work_dir = "output/life_events/test"
threads = 6

[embedding]
embedding_dim = 100
window_size = 5

[model]
model_type = "embedding"
random_state = 42
"""
    config_file.write_text(config_content)

    config = load_life_events_config(config_file)
    assert config.threads == 6
    assert config.embedding["embedding_dim"] == 100
    assert config.model["model_type"] == "embedding"


def test_load_life_events_config_env_overrides(tmp_path: Path) -> None:
    """Test LE_ prefix environment variable overrides."""
    config_file = tmp_path / "test_config.yaml"
    config_content = """
work_dir: output/life_events/test
threads: 4
"""
    config_file.write_text(config_content)

    # Set environment variables
    os.environ["LE_WORK_DIR"] = str(tmp_path / "override_work")
    os.environ["LE_THREADS"] = "8"

    try:
        config = load_life_events_config(config_file)
        assert config.work_dir == Path(tmp_path / "override_work").expanduser().resolve()
        assert config.threads == 8
    finally:
        # Clean up
        os.environ.pop("LE_WORK_DIR", None)
        os.environ.pop("LE_THREADS", None)


def test_life_events_config_defaults() -> None:
    """Test LifeEventsWorkflowConfig with default values."""
    config = LifeEventsWorkflowConfig(work_dir=Path("output/test"))
    assert config.threads == 4
    assert config.log_dir is None
    assert config.embedding == {}
    assert config.model == {}
    assert config.workflow == {}
    assert config.output == {}


def test_life_events_config_invalid_file() -> None:
    """Test error handling for missing configuration file."""
    missing_file = Path("nonexistent_config.yaml")
    with pytest.raises((FileNotFoundError, ValueError)):
        load_life_events_config(missing_file)


def test_life_events_config_malformed(tmp_path: Path) -> None:
    """Test handling of malformed YAML/JSON."""
    # Malformed YAML
    malformed_file = tmp_path / "malformed.yaml"
    malformed_file.write_text("invalid: yaml: content: [")
    with pytest.raises((ValueError, Exception)):
        load_life_events_config(malformed_file)

    # Malformed JSON
    malformed_json = tmp_path / "malformed.json"
    malformed_json.write_text("{invalid json")
    with pytest.raises((ValueError, Exception)):
        load_life_events_config(malformed_json)


def test_life_events_config_partial_sections(tmp_path: Path) -> None:
    """Test config with some sections missing."""
    config_file = tmp_path / "partial_config.yaml"
    config_content = """
work_dir: output/life_events/test
threads: 4

embedding:
  embedding_dim: 100
  window_size: 5
"""
    config_file.write_text(config_content)

    config = load_life_events_config(config_file)
    assert config.embedding["embedding_dim"] == 100
    assert config.model == {}  # Missing section should be empty dict
    assert config.workflow == {}
    assert config.output == {}


def test_life_events_config_nested_values(tmp_path: Path) -> None:
    """Test config with nested values in sections."""
    config_file = tmp_path / "nested_config.yaml"
    config_content = """
work_dir: output/life_events/test
embedding:
  embedding_dim: 100
  method: skipgram
model:
  model_type: embedding
  task_type: classification
  random_state: 42
workflow:
  save_model: true
  save_embeddings: true
output:
  format: json
  include_probabilities: true
  output_dir: output/results
"""
    config_file.write_text(config_content)

    config = load_life_events_config(config_file)
    assert config.embedding["method"] == "skipgram"
    assert config.model["random_state"] == 42
    assert config.workflow["save_embeddings"] is True
    assert config.output["output_dir"] == "output/results"
