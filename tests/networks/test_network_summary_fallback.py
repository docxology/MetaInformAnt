"""Degree statistics from the NumPy-free branch of ``get_network_summary`` match the NumPy branch.

Zero-mocks: the fallback runs in a fresh interpreter in which ``import numpy``
raises, exactly as on an install without NumPy.
"""

from __future__ import annotations

import json
import subprocess
import sys

import networkx as nx
import pytest

from metainformant.networks.analysis import graph_core

_WITHOUT_NUMPY = """
import json
import sys

sys.modules["numpy"] = None

import networkx as nx

from metainformant.networks.analysis import graph_core

assert not graph_core.HAS_NUMPY
print(json.dumps(graph_core.get_network_summary(nx.path_graph(4))["degree_stats"]))
"""


def test_degree_statistics_without_numpy_match_numpy() -> None:
    proc = subprocess.run([sys.executable, "-c", _WITHOUT_NUMPY], capture_output=True, text=True, timeout=180)
    assert proc.returncode == 0, proc.stderr[-2000:]
    fallback = json.loads(proc.stdout.strip().splitlines()[-1])

    reference = graph_core.get_network_summary(nx.path_graph(4))["degree_stats"]
    # Degrees 1, 2, 2, 1: median 1.5 and population standard deviation 0.5.
    assert reference["median"] == 1.5 and reference["std"] == pytest.approx(0.5)
    assert fallback == pytest.approx(reference)
