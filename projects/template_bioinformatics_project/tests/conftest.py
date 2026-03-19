"""
conftest.py — Shared pytest fixtures for the template bioinformatics project.

All fixtures use real file I/O on temporary directories (Zero-Mock policy).
"""

import shutil
import textwrap
from pathlib import Path

import pytest
import yaml


@pytest.fixture()
def project_root() -> Path:
    """Return the project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture()
def tmp_project(tmp_path: Path) -> Path:
    """
    Create a minimal self-contained project copy inside a temp directory.

    Copies the config/ directory and creates required subdirectories so
    scripts can be invoked without touching the real data/ or results/ trees.
    """
    project_root = Path(__file__).parent.parent

    # Config
    (tmp_path / "config").mkdir()
    config = {
        "metadata": {"project_name": "test_project", "version": "0.0.1"},
        "paths": {
            "data_raw": str(tmp_path / "data" / "raw") + "/",
            "data_processed": str(tmp_path / "data" / "processed") + "/",
            "results": str(tmp_path / "results") + "/",
            "results_figures": str(tmp_path / "results" / "figures") + "/",
            "results_tables": str(tmp_path / "results" / "tables") + "/",
            "logs": str(tmp_path / "logs") + "/",
            "config": str(tmp_path / "config" / "default.yaml"),
        },
        "processing": {
            "threads": 1,
            "filtering_threshold": 10.0,
            "normalize": True,
            "min_sample_count": 2,
            "missing_fraction_max": 0.5,
        },
        "analysis": {
            "method": "summary",
            "n_components": 3,
            "significance_threshold": 0.05,
            "correction": "bonferroni",
        },
        "visualization": {
            "dpi": 72,
            "figure_width": 6,
            "figure_height": 4,
            "color_palette": "viridis",
            "file_format": "png",
        },
        "logging": {
            "level": "DEBUG",
            "format": "%(asctime)s [%(levelname)s] %(name)s — %(message)s",
        },
    }
    config_path = tmp_path / "config" / "default.yaml"
    with config_path.open("w") as fh:
        yaml.dump(config, fh, sort_keys=False)

    # Required directories
    for d in ["data/raw", "data/processed", "results/figures", "results/tables", "logs"]:
        (tmp_path / d).mkdir(parents=True, exist_ok=True)

    return tmp_path


@pytest.fixture()
def sample_raw_csv(tmp_project: Path) -> Path:
    """Write a minimal synthetic CSV into the tmp project's data/raw/."""
    raw_dir = tmp_project / "data" / "raw"
    csv_path = raw_dir / "test_samples.csv"
    content = textwrap.dedent("""\
        sample_id,group,feature_01,feature_02,feature_03
        S0001,control,1.0,2.0,3.0
        S0002,treatment_A,1.5,2.5,3.5
        S0003,control,0.8,1.9,2.8
        S0004,treatment_B,1.2,2.2,3.2
        S0005,control,0.9,2.1,3.1
        S0006,treatment_A,1.3,2.3,3.3
        S0007,treatment_B,,2.0,3.0
        S0008,control,1.1,2.1,
    """)
    csv_path.write_text(content)
    return csv_path
