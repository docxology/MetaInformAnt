"""Optional-dependency fallbacks of the spatial modules.

Each spatial module binds its optional imports to ``None`` when the dependency is missing, so importing
the package never fails and the public functions raise ``ImportError`` with install guidance instead.
The type-checker guards around those fallbacks (``if not TYPE_CHECKING``) must not change that runtime
behaviour, so these tests run a child interpreter in which numpy, scipy, matplotlib and PIL cannot be
imported.

All tests use real implementations (real-implementation policy).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

# module -> names that must be bound to None when the optional dependencies are absent
FALLBACK_NAMES: dict[str, list[str]] = {
    "metainformant.spatial.io.merfish": ["np"],
    "metainformant.spatial.io.visium": ["np", "sp_sparse", "mmread", "Image"],
    "metainformant.spatial.io.xenium": ["np", "sp_sparse", "mmread"],
    "metainformant.spatial.analysis.neighborhood": ["np", "sp_sparse", "KDTree"],
    "metainformant.spatial.analysis.deconvolution": ["np", "sp_sparse", "nnls"],
    "metainformant.spatial.analysis.clustering": ["np", "sp_sparse", "Delaunay", "KDTree"],
    "metainformant.spatial.analysis.autocorrelation": ["np", "sp_sparse", "KDTree", "norm"],
    "metainformant.spatial.analysis.spatially_variable_genes": ["np", "false_discovery_control"],
    "metainformant.spatial.integration.scrna_mapping": ["np", "sp_sparse", "pearsonr", "spearmanr"],
    "metainformant.spatial.deconvolution.spatial_deconvolution": ["np", "scipy_nnls", "cdist"],
    "metainformant.spatial.communication.cell_communication": ["np", "cdist", "mannwhitneyu"],
    "metainformant.spatial.visualization.plots": ["np", "plt", "PatchCollection", "Wedge", "mcolors", "sp_sparse"],
}

# Modules that also expose availability flags next to the fallback names.
FLAG_MODULES = [
    "metainformant.spatial.deconvolution.spatial_deconvolution",
    "metainformant.spatial.communication.cell_communication",
]

# label -> (substring the ImportError message must contain)
EXPECTED_IMPORT_ERRORS = {
    "morans_i": "NumPy is required",
    "correlation_mapping": "NumPy is required",
    "plot_spatial_scatter": "matplotlib is required",
    "deconvolve_spots": "numpy is required",
    "compute_ligand_receptor_interactions": "numpy is required",
}

MARKER = "SPATIAL-FALLBACK-REPORT:"

CHILD_SCRIPT = """
import importlib
import json
import sys

for blocked in ("numpy", "scipy", "matplotlib", "PIL"):
    sys.modules[blocked] = None  # a None entry makes ``import <name>`` raise ImportError

names = json.loads(sys.argv[1])
flag_modules = json.loads(sys.argv[2])

bound_to_none = {}
for module_name, attrs in names.items():
    module = importlib.import_module(module_name)
    bound_to_none[module_name] = {attr: getattr(module, attr) is None for attr in attrs}

flags = {}
for module_name in flag_modules:
    module = importlib.import_module(module_name)
    flags[module_name] = {flag: getattr(module, flag) for flag in ("HAS_NUMPY", "HAS_SCIPY")}

from metainformant.spatial.analysis.autocorrelation import morans_i
from metainformant.spatial.communication.cell_communication import compute_ligand_receptor_interactions
from metainformant.spatial.deconvolution.spatial_deconvolution import deconvolve_spots
from metainformant.spatial.integration.scrna_mapping import correlation_mapping
from metainformant.spatial.visualization.plots import plot_spatial_scatter

calls = {
    "morans_i": lambda: morans_i([1.0, 2.0, 3.0], None),
    "correlation_mapping": lambda: correlation_mapping([[1.0, 2.0]], [[1.0, 2.0]]),
    "plot_spatial_scatter": lambda: plot_spatial_scatter([[0.0, 0.0]], [1.0], "unused.png"),
    "deconvolve_spots": lambda: deconvolve_spots([[1.0]], {"A": [1.0]}),
    "compute_ligand_receptor_interactions": lambda: compute_ligand_receptor_interactions([[1.0]], ["A"]),
}
import_errors = {}
for label, call in calls.items():
    try:
        call()
    except ImportError as exc:
        import_errors[label] = str(exc)
    else:
        import_errors[label] = None

print(sys.argv[3] + json.dumps({"bound_to_none": bound_to_none, "flags": flags, "import_errors": import_errors}))
"""


@pytest.fixture(scope="module")
def blocked_report(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """Import every spatial module in a child interpreter that cannot import the numeric stack."""
    workdir: Path = tmp_path_factory.mktemp("blocked_imports")
    proc = subprocess.run(
        [sys.executable, "-c", CHILD_SCRIPT, json.dumps(FALLBACK_NAMES), json.dumps(FLAG_MODULES), MARKER],
        capture_output=True,
        text=True,
        timeout=180,
        cwd=workdir,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    line = next(line for line in proc.stdout.splitlines() if line.startswith(MARKER))
    report: dict[str, Any] = json.loads(line[len(MARKER) :])
    return report


class TestMissingDependencies:
    def test_fallback_names_are_bound_to_none(self, blocked_report: dict[str, Any]) -> None:
        unbound = {
            module: [name for name, is_none in names.items() if not is_none]
            for module, names in blocked_report["bound_to_none"].items()
        }
        assert {module: names for module, names in unbound.items() if names} == {}

    def test_availability_flags_are_false(self, blocked_report: dict[str, Any]) -> None:
        for module in FLAG_MODULES:
            assert blocked_report["flags"][module] == {"HAS_NUMPY": False, "HAS_SCIPY": False}, module

    @pytest.mark.parametrize("label", sorted(EXPECTED_IMPORT_ERRORS))
    def test_public_functions_raise_import_error_with_guidance(
        self, blocked_report: dict[str, Any], label: str
    ) -> None:
        message = blocked_report["import_errors"][label]
        assert message is not None, f"{label} did not raise ImportError"
        assert EXPECTED_IMPORT_ERRORS[label] in message
        assert "uv pip install" in message
