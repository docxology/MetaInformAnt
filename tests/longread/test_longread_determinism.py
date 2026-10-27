"""Reproducibility checks for the long-read package.

Audit result: nothing under ``metainformant.longread`` draws random numbers, so no analysis path is
unseeded. The sources of run-to-run variation that were found are hash-order dependent set
iteration (fixed in ``analysis.structural`` and ``utils.summary``, tested in
``test_longread_structural_merge.py``) and the completion-order keys of parallel batch results
(tested here). The AST guard below keeps the first finding true: random number generation has to be
added deliberately, with a seed, rather than slipping in unseeded.
"""

from __future__ import annotations

import ast
import time
from pathlib import Path

import metainformant.longread
from metainformant.longread.utils.batch import batch_compute_metrics, process_batch

PACKAGE_DIR = Path(metainformant.longread.__file__).resolve().parent

RANDOM_MODULES = {"random", "secrets", "uuid"}
RANDOM_NAMES = {"default_rng", "RandomState", "urandom", "SeedSequence", "shuffle", "randint", "uniform"}


def _random_usages(source: str) -> list[str]:
    """Names of random-number facilities a module imports or uses, found from its syntax tree."""
    found: list[str] = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            found += [alias.name for alias in node.names if alias.name.split(".")[0] in RANDOM_MODULES]
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module.split(".")[0] in RANDOM_MODULES or node.module.endswith("random"):
                found.append(node.module)
        elif isinstance(node, ast.Attribute) and (node.attr == "random" or node.attr in RANDOM_NAMES):
            found.append(node.attr)
    return found


class TestNoUnseededRandomness:
    def test_the_detector_finds_random_usage(self) -> None:
        # the guard itself works on the usual spellings
        assert _random_usages("import random\nrandom.random()")
        assert _random_usages("import numpy as np\nnp.random.rand(3)")
        assert _random_usages("import numpy as np\nrng = np.random.default_rng()")
        assert _random_usages("from numpy.random import default_rng")
        assert _random_usages("import secrets")
        assert not _random_usages("import math\nx = math.sqrt(2)")

    def test_no_module_of_the_package_uses_random_number_generation(self) -> None:
        offenders = {}
        for path in sorted(PACKAGE_DIR.rglob("*.py")):
            usages = _random_usages(path.read_text(encoding="utf-8"))
            if usages:
                offenders[str(path.relative_to(PACKAGE_DIR))] = usages
        assert not offenders, (
            "random number generation found; every analysis path must take an explicit seed or be deterministic: "
            f"{offenders}"
        )


class TestParallelBatchOrder:
    ITEMS = {f"item{i}": i for i in range(6)}

    @classmethod
    def _slow_in_input_order(cls, item_id: str, value: int) -> int:
        time.sleep(0.04 * (len(cls.ITEMS) - value))  # later items finish first
        if value in (1, 4):
            raise ValueError(f"cannot process {item_id}")
        return value * 10

    def test_results_and_errors_follow_the_input_order_not_the_completion_order(self) -> None:
        result = process_batch(self.ITEMS, self._slow_in_input_order, max_workers=4)
        assert list(result.results) == ["item0", "item2", "item3", "item5"]
        assert list(result.errors) == ["item1", "item4"]
        assert result.results == {"item0": 0, "item2": 20, "item3": 30, "item5": 50}
        assert (result.successful, result.failed, result.total_items) == (4, 2, 6)

    def test_sequential_and_parallel_runs_agree(self) -> None:
        sequential = process_batch(self.ITEMS, self._slow_in_input_order)
        parallel = process_batch(self.ITEMS, self._slow_in_input_order, max_workers=3)
        assert list(sequential.results) == list(parallel.results)
        assert list(sequential.errors) == list(parallel.errors)
        assert sequential.results == parallel.results

    def test_batch_metrics_report_samples_in_input_order(self) -> None:
        read_sets = {f"sample{i}": [1000 * (i + 1)] * (6 - i) for i in range(6)}
        result = batch_compute_metrics(read_sets, max_workers=4)
        assert list(result.results) == list(read_sets)
