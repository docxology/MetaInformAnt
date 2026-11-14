"""Cold producer startup loads the streaming producer without loading plotting."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_producer_cold_start_does_not_initialize_plotting(tmp_path: Path) -> None:
    """Import the producer through the real import chain with a fresh cache within 25s."""
    entry = "metainformant.rna.engine.streaming_orchestrator"
    env = dict(os.environ)
    for key, name in [("HOME", "home"), ("TMPDIR", "tmp"), ("MPLCONFIGDIR", "mpl"), ("XDG_CACHE_HOME", "cache")]:
        path = tmp_path / name
        path.mkdir()
        env[key] = str(path)
    env["PYTHONPATH"] = str(ROOT / "src")
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
    """Every exported module is the same object through attribute access and import."""
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
