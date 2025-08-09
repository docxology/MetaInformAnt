import sys
from pathlib import Path


def test_cli_run_config_smoke(monkeypatch, tmp_path: Path):
    # Create minimal config file
    cfg_text = (
        "work_dir: " + str(tmp_path / "work") + "\n"
        "threads: 1\n"
        "species_list: []\n"
        "steps: {}\n"
    )
    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text(cfg_text, encoding="utf-8")

    # Bypass actual execution by monkeypatching execute_workflow
    import metainformant.rna.workflow as wf
    monkeypatch.setattr(wf, "execute_workflow", lambda cfg, check=False: [0])

    # Run CLI main
    from metainformant.__main__ import main
    sys.argv = ["metainformant", "rna", "run-config", "--config", str(cfg_file)]
    main()


