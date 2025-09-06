from pathlib import Path

from metainformant.rna.workflow import load_workflow_config, plan_workflow


def test_load_workflow_config_yaml(tmp_path: Path):
    cfg_text = (
        "work_dir: " + str(tmp_path / "work") + "\n"
        "log_dir: " + str(tmp_path / "logs") + "\n"
        "threads: 2\n"
        "species_list:\n  - X\n"
        "steps:\n  metadata: { a: 1 }\n  quant: { threads: 2 }\n"
    )
    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text(cfg_text, encoding="utf-8")

    cfg = load_workflow_config(cfg_file)
    assert cfg.work_dir == (tmp_path / "work").resolve()
    assert cfg.log_dir == (tmp_path / "logs").resolve()
    assert cfg.threads == 2
    assert cfg.species_list == ["X"]
    steps = dict(plan_workflow(cfg))
    # common params merged with per-step
    assert "metadata" in steps and isinstance(steps["metadata"], dict)
    assert "quant" in steps and steps["quant"]["threads"] == 2
