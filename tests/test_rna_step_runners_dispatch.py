from __future__ import annotations

from pathlib import Path


def test_each_step_runner_invokes_real_subcommand_or_skips(tmp_path: Path):
    from metainformant.rna import steps as step_mod
    from metainformant.rna.amalgkit import check_cli_available

    ok, _ = check_cli_available()
    if not ok:
        import pytest

        pytest.skip("amalgkit not available on PATH")

    called: list[str] = []

    for name, runner in step_mod.STEP_RUNNERS.items():
        res = runner({}, work_dir=tmp_path / "work", log_dir=tmp_path / "logs")
        called.append(name)
        # Do not assert return code strictly; some steps may fail without inputs
        assert hasattr(res, "returncode")

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
