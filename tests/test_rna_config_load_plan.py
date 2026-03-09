"""Tests for RNA workflow configuration loading and planning.

Tests YAML config parsing, environment variable overrides, and workflow planning.
"""

from __future__ import annotations

from pathlib import Path


def test_load_workflow_config_and_plan_uses_yaml_values():
    """Test that workflow config loading and planning uses YAML values correctly."""
    from metainformant.rna.engine.workflow import load_workflow_config, plan_workflow

    repo_root = Path(__file__).resolve().parents[1]
    # Use test config file with threads: 6 instead of production config with threads: 12
    cfg_path = repo_root / "config" / "amalgkit" / "amalgkit_test.yaml"
    if not cfg_path.exists():
        # Fallback to pbarbatus config if test config doesn't exist
        cfg_path = repo_root / "config" / "amalgkit" / "amalgkit_pogonomyrmex_barbatus.yaml"
        # Adjust expectation for pbarbatus config which has threads: 12
        expected_threads = 12
    else:
        expected_threads = 6

    # Sanity: file exists in repo (either test or pbarbatus)
    assert cfg_path.exists(), f"Config file not found: {cfg_path}"

    cfg = load_workflow_config(cfg_path)

    # Core fields from YAML
    if "test" in cfg_path.name:
        assert cfg.work_dir.as_posix().endswith("output/amalgkit/test/work")
        assert (cfg.log_dir is not None) and cfg.log_dir.as_posix().endswith("output/amalgkit/test/logs")
        expected_fastq_dir = "output/amalgkit/test/fastq"
        expected_merge_out = "output/amalgkit/test/merged/merged_abundance.tsv"
    else:
        assert cfg.work_dir.as_posix().endswith("output/amalgkit/pbarbatus/work")
        assert (cfg.log_dir is not None) and cfg.log_dir.as_posix().endswith("output/amalgkit/pbarbatus/logs")
        expected_fastq_dir = "output/amalgkit/pbarbatus/fastq"
        expected_merge_out = "output/amalgkit/pbarbatus/merged/merged_abundance.tsv"

    assert cfg.threads == expected_threads

    steps = plan_workflow(cfg)
    names = [n for n, _ in steps]
    assert names[0] == "metadata" and names[-1] == "sanity"

    params = {n: p for n, p in steps}

    # Per-step params merged with common
    assert params["getfastq"].get("out_dir", "").endswith(expected_fastq_dir)
    assert params["merge"].get("out", "").endswith(expected_merge_out)
    # Common threads propagated into steps that support them (getfastq, integrate, quant)
    steps_with_threads = {"getfastq", "integrate", "quant"}
    assert all(params[n].get("threads") == expected_threads for n in steps_with_threads if n in params)


def test_env_overrides_for_config_threads(tmp_path: Path):
    """Test that environment variables can override config values."""
    from metainformant.rna.engine.workflow import load_workflow_config

    repo_root = Path(__file__).resolve().parents[1]
    # Use test config if available, otherwise use pbarbatus
    cfg_path = repo_root / "config" / "amalgkit" / "amalgkit_test.yaml"
    if not cfg_path.exists():
        cfg_path = repo_root / "config" / "amalgkit" / "amalgkit_pogonomyrmex_barbatus.yaml"

    # Ensure at least one config exists
    assert cfg_path.exists(), f"Config file not found: {cfg_path}"

    import os

    old = os.environ.get("AK_THREADS")
    try:
        os.environ["AK_THREADS"] = "3"
        cfg = load_workflow_config(cfg_path)
        assert cfg.threads == 3, f"Expected threads=3, got {cfg.threads}"
    finally:
        # Restore original value or remove if not set
        if old is None:
            os.environ.pop("AK_THREADS", None)
        else:
            os.environ["AK_THREADS"] = old
