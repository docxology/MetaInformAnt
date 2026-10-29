"""Cold worker startup uses acquisition code without loading plotting."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from tests.rna.test_ehd_runtime_paths import SOURCE_OPS, require_operator_files

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(
    "entry",
    [pytest.param("campaign", marks=pytest.mark.external_tool), "metainformant.rna.engine.acquisition_worker"],
)
def test_worker_cold_start_does_not_initialize_plotting(tmp_path: Path, entry: str) -> None:
    """Use the real import chain and a fresh cache under the original 25s budget."""
    if entry == "campaign":
        require_operator_files("campaign.py")
    env = dict(os.environ)
    for key, name in [("HOME", "home"), ("TMPDIR", "tmp"), ("MPLCONFIGDIR", "mpl"), ("XDG_CACHE_HOME", "cache")]:
        path = tmp_path / name
        path.mkdir()
        env[key] = str(path)
    env["PYTHONPATH"] = os.pathsep.join([str(ROOT / "src"), str(SOURCE_OPS)])
    code = """
import importlib,json,sys
importlib.import_module(sys.argv[1])
assert "matplotlib.font_manager" not in sys.modules
assert "matplotlib.pyplot" not in sys.modules
assert "metainformant.rna.analysis.atlas_plots" not in sys.modules
print(json.dumps({"entry":sys.argv[1],"plotting_loaded":False}))
"""
    result = subprocess.run(
        [sys.executable, "-B", "-c", code, entry], env=env, capture_output=True, text=True, timeout=25
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout) == {"entry": entry, "plotting_loaded": False}
    assert not list((tmp_path / "mpl").iterdir())


def test_rna_lazy_public_modules_preserve_identity() -> None:
    """Every existing exported module remains importable through either API."""
    import importlib

    import metainformant.rna as rna

    expected = {"amalgkit", "analysis", "core", "deconvolution", "engine", "retrieval", "splicing", "steps"}
    assert set(rna.__all__) == expected
    assert expected.issubset(dir(rna))
    for name in sorted(expected):
        module = getattr(rna, name)
        assert module is importlib.import_module(f"metainformant.rna.{name}")
        assert getattr(rna, name) is module


def test_rna_unknown_public_attribute_refuses() -> None:
    import metainformant.rna as rna

    with pytest.raises(AttributeError, match="does_not_exist"):
        getattr(rna, "does_not_exist")
