from __future__ import annotations

from pathlib import Path


def test_load_workflow_config_and_plan_uses_yaml_values():
    from metainformant.rna.workflow import load_workflow_config, plan_workflow

    repo_root = Path(__file__).resolve().parents[1]
    cfg_path = repo_root / "config" / "amalgkit_pbarbatus.yaml"

    # Sanity: file exists in repo
    assert cfg_path.exists()

    cfg = load_workflow_config(cfg_path)

    # Core fields from YAML
    assert cfg.work_dir.as_posix().endswith("output/amalgkit/pbarbatus/work")
    assert (cfg.log_dir is not None) and cfg.log_dir.as_posix().endswith("output/amalgkit/pbarbatus/logs")
    assert cfg.threads == 8

    steps = plan_workflow(cfg)
    names = [n for n, _ in steps]
    assert names[0] == "metadata" and names[-1] == "sanity"

    params = {n: p for n, p in steps}

    # Per-step params merged with common
    assert params["getfastq"].get("out-dir", "").endswith("output/amalgkit/pbarbatus/fastq")
    assert params["merge"].get("out", "").endswith("output/amalgkit/pbarbatus/merged/merged_abundance.tsv")
    # Common threads propagated into each step
    assert all(p.get("threads") == 8 for p in params.values())


def test_env_overrides_for_config_threads(tmp_path: Path):
    from metainformant.rna.workflow import load_workflow_config

    repo_root = Path(__file__).resolve().parents[1]
    cfg_path = repo_root / "config" / "amalgkit_pbarbatus.yaml"

    import os
    old = os.environ.get("AK_THREADS")
    os.environ["AK_THREADS"] = "3"
    cfg = load_workflow_config(cfg_path)
    assert cfg.threads == 3
    if old is None:
        del os.environ["AK_THREADS"]
    else:
        os.environ["AK_THREADS"] = old


