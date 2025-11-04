"""Life events workflow configuration management."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..core.config import apply_env_overrides, load_mapping_from_file


@dataclass
class LifeEventsWorkflowConfig:
    """Configuration for life events workflow execution."""

    work_dir: Path
    threads: int = 4
    log_dir: Path | None = None
    embedding: dict[str, Any] = field(default_factory=dict)
    model: dict[str, Any] = field(default_factory=dict)
    workflow: dict[str, Any] = field(default_factory=dict)
    output: dict[str, Any] = field(default_factory=dict)


def load_life_events_config(config_file: str | Path) -> LifeEventsWorkflowConfig:
    """Load LifeEventsWorkflowConfig from a config file with env overrides.

    Expected top-level keys:
      - work_dir (str): Working directory for outputs
      - log_dir (str, optional): Directory for logs
      - threads (int, default: 4): Number of threads for parallel processing
      - embedding (mapping, optional): Embedding configuration
        - embedding_dim (int, default: 100)
        - window_size (int, default: 5)
        - epochs (int, default: 10)
        - method (str, default: "skipgram")
        - learning_rate (float, default: 0.01)
      - model (mapping, optional): Model configuration
        - model_type (str, default: "embedding")
        - task_type (str, default: "classification")
        - random_state (int, optional)
      - workflow (mapping, optional): Workflow-specific settings
      - output (mapping, optional): Output configuration

    Environment variables can override config values using prefix "LE":
      - LE_WORK_DIR
      - LE_THREADS
      - LE_EMBEDDING_DIM
      - LE_WINDOW_SIZE
      - etc.

    Args:
        config_file: Path to configuration file (YAML, TOML, or JSON)

    Returns:
        LifeEventsWorkflowConfig instance

    Examples:
        >>> config = load_life_events_config("config/life_events.yaml")
        >>> print(config.work_dir)
        >>> print(config.embedding["embedding_dim"])
    """
    raw = load_mapping_from_file(config_file)
    raw = apply_env_overrides(raw, prefix="LE")

    work_dir = Path(raw.get("work_dir", "output/life_events/work")).expanduser().resolve()
    log_dir_val = raw.get("log_dir")
    log_dir = Path(log_dir_val).expanduser().resolve() if isinstance(log_dir_val, str) else None
    threads = int(raw.get("threads", 4))

    embedding_cfg = raw.get("embedding") if isinstance(raw.get("embedding"), dict) else {}
    model_cfg = raw.get("model") if isinstance(raw.get("model"), dict) else {}
    workflow_cfg = raw.get("workflow") if isinstance(raw.get("workflow"), dict) else {}
    output_cfg = raw.get("output") if isinstance(raw.get("output"), dict) else {}

    return LifeEventsWorkflowConfig(
        work_dir=work_dir,
        log_dir=log_dir,
        threads=threads,
        embedding=embedding_cfg,
        model=model_cfg,
        workflow=workflow_cfg,
        output=output_cfg,
    )

