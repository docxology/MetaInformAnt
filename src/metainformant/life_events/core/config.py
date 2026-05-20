"""Configuration management for life events module.

This module provides configuration classes and loading functions for life events workflows.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional

from metainformant.core.utils import config as core_config
from metainformant.core.utils import logging


@dataclass
class LifeEventsWorkflowConfig:
    """Configuration for life events workflow processing.

    Attributes:
        work_dir: Working directory for outputs
        embedding_dim: Dimensionality for event embeddings
        model_type: Type of prediction model ('embedding', 'lstm', 'ensemble')
        learning_rate: Learning rate for model training
        batch_size: Batch size for training
        epochs: Number of training epochs
        random_seed: Random seed for reproducibility
        output_dir: Directory for output files
        log_level: Logging level
        use_gpu: Whether to use GPU if available
        save_models: Whether to save trained models
        cross_validation_folds: Number of CV folds
        test_split: Fraction of data for testing
        embedding_window: Context window size for embeddings
        min_event_count: Minimum event frequency for embedding
        sequence_max_length: Maximum sequence length for models
        threads: Number of threads for parallel processing
        embedding: Nested embedding section from config files
        model: Nested model section from config files
        workflow: Nested workflow section from config files
        output: Nested output section from config files
    """

    work_dir: Path
    embedding_dim: int = 100
    model_type: str = "embedding"
    learning_rate: float = 0.001
    batch_size: int = 32
    epochs: int = 100
    random_seed: int = 42
    output_dir: Optional[Path] = None
    log_level: str = "INFO"
    use_gpu: bool = False
    save_models: bool = True
    cross_validation_folds: int = 5
    test_split: float = 0.2
    embedding_window: int = 5
    min_event_count: int = 5
    sequence_max_length: int = 100
    threads: int = 4
    embedding: Dict[str, Any] = field(default_factory=dict)
    model: Dict[str, Any] = field(default_factory=dict)
    workflow: Dict[str, Any] = field(default_factory=dict)
    output: Dict[str, Any] = field(default_factory=dict)
    log_dir: Optional[Path] = None

    def __post_init__(self):
        """Post-initialization validation and setup."""
        self.work_dir = self.work_dir.expanduser().resolve()

        # Older callers passed scalar aliases. Preserve them while normalizing
        # the public config sections to dictionaries.
        if not isinstance(self.embedding, dict):
            self.embedding_dim = int(self.embedding)
            self.embedding = {"embedding_dim": self.embedding_dim}
        if not isinstance(self.model, dict):
            self.model_type = str(self.model)
            self.model = {"model_type": self.model_type}
        if not isinstance(self.workflow, dict):
            self.workflow = {}
        if not isinstance(self.output, dict):
            self.output = {}

        if self.embedding:
            self.embedding_dim = int(self.embedding.get("embedding_dim", self.embedding_dim))
            self.embedding_window = int(self.embedding.get("window_size", self.embedding_window))
            self.epochs = int(self.embedding.get("epochs", self.epochs))
            self.learning_rate = float(self.embedding.get("learning_rate", self.learning_rate))
        if self.model:
            self.model_type = str(self.model.get("model_type", self.model_type))
            if "random_state" in self.model:
                self.random_seed = int(self.model["random_state"])

        if self.log_dir is not None and self.output_dir is None:
            self.log_dir = self.log_dir.expanduser().resolve()
            self.output_dir = self.log_dir

        if self.output_dir is None:
            output_dir = self.output.get("output_dir")
            self.output_dir = Path(output_dir).expanduser().resolve() if output_dir else self.work_dir / "output"
        else:
            self.output_dir = self.output_dir.expanduser().resolve()

        # Ensure directories exist
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary."""
        return {
            "work_dir": str(self.work_dir),
            "embedding_dim": self.embedding_dim,
            "model_type": self.model_type,
            "learning_rate": self.learning_rate,
            "batch_size": self.batch_size,
            "epochs": self.epochs,
            "random_seed": self.random_seed,
            "output_dir": str(self.output_dir) if self.output_dir else None,
            "log_level": self.log_level,
            "use_gpu": self.use_gpu,
            "save_models": self.save_models,
            "cross_validation_folds": self.cross_validation_folds,
            "test_split": self.test_split,
            "embedding_window": self.embedding_window,
            "min_event_count": self.min_event_count,
            "sequence_max_length": self.sequence_max_length,
            "threads": self.threads,
            "embedding": self.embedding,
            "model": self.model,
            "workflow": self.workflow,
            "output": self.output,
            "log_dir": str(self.log_dir) if self.log_dir else None,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> LifeEventsWorkflowConfig:
        """Create config from dictionary."""
        data = dict(data)

        # Convert string paths to Path objects
        if "work_dir" in data:
            data["work_dir"] = Path(data["work_dir"]).expanduser().resolve()
        if "output_dir" in data and data["output_dir"]:
            data["output_dir"] = Path(data["output_dir"]).expanduser().resolve()
        if "log_dir" in data and data["log_dir"]:
            data["log_dir"] = Path(data["log_dir"]).expanduser().resolve()

        if "threads" in data:
            try:
                data["threads"] = int(data["threads"])
            except (TypeError, ValueError):
                data["threads"] = cls.__dataclass_fields__["threads"].default

        return cls(**data)


def load_life_events_config(config_file: str | Path, prefix: str = "LE") -> LifeEventsWorkflowConfig:
    """Load life events configuration from file with environment overrides.

    Args:
        config_file: Path to configuration file (YAML/TOML/JSON)
        prefix: Environment variable prefix for overrides

    Returns:
        LifeEventsWorkflowConfig instance

    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If config is invalid
    """
    # Load raw config
    raw_config = core_config.load_mapping_from_file(config_file)

    # Apply environment variable overrides
    raw_config = core_config.apply_env_overrides(raw_config, prefix=prefix)

    # Validate required fields
    if "work_dir" not in raw_config:
        raise ValueError("Configuration must specify 'work_dir'")

    # Convert and validate
    try:
        config = LifeEventsWorkflowConfig.from_dict(raw_config)
    except Exception as e:
        raise ValueError(f"Invalid configuration: {e}")

    return config


def create_default_config(work_dir: str | Path) -> LifeEventsWorkflowConfig:
    """Create a default configuration for life events workflows.

    Args:
        work_dir: Working directory path

    Returns:
        Default LifeEventsWorkflowConfig
    """
    return LifeEventsWorkflowConfig(work_dir=Path(work_dir))


def save_config(config: LifeEventsWorkflowConfig, output_file: str | Path) -> None:
    """Save configuration to file.

    Args:
        config: Configuration to save
        output_file: Output file path
    """
    config_dict = config.to_dict()

    # Determine format from extension
    if str(output_file).endswith(".yaml") or str(output_file).endswith(".yml"):
        import yaml

        with open(output_file, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False)
    elif str(output_file).endswith(".json"):
        import json

        with open(output_file, "w") as f:
            json.dump(config_dict, f, indent=2)
    else:
        raise ValueError(f"Unsupported config format: {output_file}")

    logger = logging.get_logger(__name__)
    logger.info(f"Saved configuration to {output_file}")
