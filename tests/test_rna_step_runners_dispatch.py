from __future__ import annotations

from pathlib import Path


def test_each_step_runner_invokes_correct_subcommand(monkeypatch, tmp_path: Path):
    from metainformant.rna import steps as step_mod
    import metainformant.rna.amalgkit as ak

    called: list[str] = []

    def fake_run_amalgkit(subcommand, params=None, **kwargs):
        called.append(subcommand)
        class R:
            returncode = 0
            stdout = ""
            stderr = ""
        return R()

    monkeypatch.setattr(ak, "run_amalgkit", fake_run_amalgkit)

    for name, runner in step_mod.STEP_RUNNERS.items():
        res = runner({}, work_dir=tmp_path / "work", log_dir=tmp_path / "logs")
        assert res.returncode == 0

    # All expected names were called at least once
    expected = {
        "metadata",
        "integrate",
        "config",
        "select",
        "getfastq",
        "quant",
        "merge",
        "cstmm",
        "curate",
        "csca",
        "sanity",
    }
    assert expected.issubset(set(called))


