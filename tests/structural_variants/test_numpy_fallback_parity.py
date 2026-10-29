"""Regression tests for the optional-numpy branches of the structural-variant modules.

``merge``, ``sv_calling``, ``cnv`` and ``plots`` treat numpy as an optional
dependency and carry a pure-Python branch for installs without it. Both branches
must compute the same statistics. The numpy branch uses the true median (the mean
of the two middle values for an even count), a sample standard deviation that
defaults to 100 for a single observation, and median-centred log2 ratios with
optional GC correction. The pure-Python branches used to take the upper median,
reported a zero standard deviation for a single insert size and returned raw,
uncorrected log2 ratios.

Every scenario runs the real code twice: in this process (numpy available) and in
a fresh interpreter in which the four modules are imported while ``import numpy``
raises ``ImportError``, so each module's optional import falls back to ``None``.
matplotlib is loaded first (it needs numpy to draw) and numpy is made importable
again once the modules are loaded; only the modules' own optional imports fail.
No test doubles are involved.
"""

from __future__ import annotations

import json
import os
import random
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

_RESULT_MARKER = "NUMPY_FALLBACK_RESULTS="

# Position groups with odd and even sizes; the true medians are 150, 2.5 -> 2,
# 5, 2 and 25 (the integer value is truncated like ``int(np.median(...))``).
_POSITION_GROUPS = ([100, 200], [1, 4], [5], [3, 1, 2], [10, 20, 30, 41])
_EXPECTED_POSITION_MEDIANS = [150, 2, 5, 2, 25]

# Insert-size samples: a single observation, two, four and six.
_INSERT_SIZE_SAMPLES = ([400], [300, 500], [300, 400, 500, 700], [250, 260, 270, 900, 950, 111])
_EXPECTED_INSERT_SIZE_STATS = [
    # mean, std, median, mad
    [400.0, 100.0, 400.0, 0.0],  # single observation: default std of 100
    [400.0, 141.4213562373095, 400.0, 100.0],
    [475.0, 170.78251276599332, 450.0, 100.0],
    [456.8333333333333, 367.5978327828752, 265.0, 84.5],
]

# Four discordant pairs (inward-facing, ~4 kb apart) forming one deletion cluster.
_DISCORDANT_PAIRS = [(1000, 5000), (1010, 4990), (1020, 4995), (1031, 5001)]


def collect_results() -> dict[str, Any]:
    """Run the numpy-optional code paths and return JSON-serialisable results.

    Imports happen inside the function so that a subprocess can hide numpy first.
    """
    from metainformant.structural_variants.detection import cnv, sv_calling
    from metainformant.structural_variants.filtering import merge
    from metainformant.structural_variants.visualization import plots

    results: dict[str, Any] = {"numpy_available": [m.np is not None for m in (merge, sv_calling, cnv, plots)]}

    results["median_position"] = [merge._median_position(list(group)) for group in _POSITION_GROUPS]

    merged, _ = merge.merge_callsets(
        {
            "caller_a": [{"chrom": "chr1", "start": 10000, "end": 15000, "sv_type": "DEL"}],
            "caller_b": [{"chrom": "chr1", "start": 10100, "end": 15200, "sv_type": "DEL"}],
        },
        min_overlap=0.5,
    )
    results["merge_callsets"] = [[m.chrom, m.start, m.end, m.sv_type, m.size] for m in merged]

    # The same translocation written from opposite sides by two callers (second breakpoints 20000 and 20101).
    junction_merged, _ = merge.merge_callsets(
        {
            "caller_a": [
                {"chrom": "chr1", "start": 10000, "end": 10000, "sv_type": "TRA", "chrom2": "chr2", "pos2": 20000}
            ],
            "caller_b": [
                {"chrom": "chr2", "start": 20101, "end": 20101, "sv_type": "TRA", "chrom2": "chr1", "pos2": 10100}
            ],
        },
        max_distance=500,
    )
    results["merge_callsets_translocation"] = [[m.chrom, m.start, m.chrom2, m.pos2] for m in junction_merged]

    alignments = [
        {
            "name": f"pair{i}",
            "chrom": "chr1",
            "pos": pos,
            "mapq": 60,
            "cigar": "100M",
            "mate_chrom": "chr1",
            "mate_pos": mate_pos,
            "insert_size": mate_pos - pos + 100,
            "is_reverse": False,
            "mate_is_reverse": True,
        }
        for i, (pos, mate_pos) in enumerate(_DISCORDANT_PAIRS)
    ]
    calls = sv_calling.call_structural_variants(
        alignments,
        min_support=3,
        insert_size_stats=sv_calling.InsertSizeStats(mean=400.0, std=50.0),
    )
    results["sv_calls"] = [[c.chrom, c.start, c.end, c.sv_type.value, c.size] for c in calls]

    insert_stats = []
    for sizes in _INSERT_SIZE_SAMPLES:
        reads = [{"chrom": "chr1", "mate_chrom": "chr1", "insert_size": size} for size in sizes]
        estimate = sv_calling._estimate_insert_size(reads)
        insert_stats.append([estimate.mean, estimate.std, estimate.median, estimate.mad])
    results["insert_size_stats"] = insert_stats

    rng = random.Random(5)
    tumor = [rng.randint(5, 80) for _ in range(60)]
    normal = [rng.randint(5, 60) for _ in range(60)]
    gc = [rng.random() for _ in range(60)]
    gc[3], gc[4], gc[5], gc[6] = 0.15, 0.45, 1.0, 0.0  # values on the GC bin edges
    results["log2_ratio"] = [float(v) for v in cnv.calculate_log2_ratio(tumor, normal)]
    results["log2_ratio_gc"] = [float(v) for v in cnv.calculate_log2_ratio(tumor, normal, gc_content=gc)]

    figure = plots.plot_coverage_track([10.0, 20.0, 0.0, 0.0], [], {"chrom": "chr1", "start": 0, "end": 4000})
    results["coverage_median_label"] = figure.axes[0].get_legend().get_texts()[0].get_text()

    return results


# Load matplotlib and its Agg backend first (they need numpy), import the modules under test while
# numpy is hidden so their optional import falls back to None, then make numpy importable again.
_RUNNER = (
    "import runpy, sys\n"
    "import matplotlib\n"
    "matplotlib.use('Agg')\n"
    "import matplotlib.pyplot as plt\n"
    "plt.close(plt.figure())\n"
    "real_numpy = sys.modules['numpy']\n"
    "sys.modules['numpy'] = None\n"
    "import metainformant.structural_variants.detection.cnv\n"
    "import metainformant.structural_variants.detection.sv_calling\n"
    "import metainformant.structural_variants.filtering.merge\n"
    "import metainformant.structural_variants.visualization.plots\n"
    "sys.modules['numpy'] = real_numpy\n"
    "runpy.run_path(sys.argv[1], run_name='__main__')\n"
)


def _run_without_numpy(workdir: Path) -> dict[str, Any]:
    """Execute :func:`collect_results` in a fresh interpreter whose modules were imported without numpy."""
    import metainformant.structural_variants as structural_variants

    # The child must import the very copy of the package that this process is testing: the directory that holds
    # ``metainformant`` (``src`` in the repository), not whichever copy the environment lists first.
    package_root = Path(structural_variants.__file__).parents[2]
    env = dict(os.environ, PYTHONPATH=str(package_root))
    env["MPLBACKEND"] = "Agg"
    completed = subprocess.run(
        [sys.executable, "-c", _RUNNER, str(Path(__file__).resolve())],
        capture_output=True,
        text=True,
        cwd=workdir,
        env=env,
        timeout=300,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr[-4000:]
    lines = [line for line in completed.stdout.splitlines() if line.startswith(_RESULT_MARKER)]
    assert len(lines) == 1, completed.stdout[-2000:]
    parsed: dict[str, Any] = json.loads(lines[0][len(_RESULT_MARKER) :])
    return parsed


@pytest.fixture(scope="module")
def numpy_results() -> dict[str, Any]:
    return collect_results()


@pytest.fixture(scope="module")
def fallback_results(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    return _run_without_numpy(tmp_path_factory.mktemp("no_numpy"))


def _assert_rows_close(actual: list[list[float]], expected: list[list[float]]) -> None:
    assert len(actual) == len(expected)
    for actual_row, expected_row in zip(actual, expected):
        assert actual_row == pytest.approx(expected_row, rel=1e-9, abs=1e-9)


def test_fallback_interpreter_really_lacks_numpy(
    numpy_results: dict[str, Any], fallback_results: dict[str, Any]
) -> None:
    assert numpy_results["numpy_available"] == [True, True, True, True]
    assert fallback_results["numpy_available"] == [False, False, False, False]


def test_merge_median_position_is_the_true_median(
    numpy_results: dict[str, Any], fallback_results: dict[str, Any]
) -> None:
    assert numpy_results["median_position"] == _EXPECTED_POSITION_MEDIANS
    assert fallback_results["median_position"] == _EXPECTED_POSITION_MEDIANS


def test_merge_callsets_consensus_is_the_median_of_both_callers(
    numpy_results: dict[str, Any], fallback_results: dict[str, Any]
) -> None:
    # starts 10000/10100 -> 10050, ends 15000/15200 -> 15100
    expected = [["chr1", 10050, 15100, "DEL", 5050]]
    assert numpy_results["merge_callsets"] == expected
    assert fallback_results["merge_callsets"] == expected


def test_merged_translocation_breakpoints_are_the_median_of_each_breakpoint(
    numpy_results: dict[str, Any], fallback_results: dict[str, Any]
) -> None:
    # first breakpoints 10000 and 10100 -> 10050; second breakpoints 20000 and 20101 -> 20050.5 -> 20050 (truncated,
    # as for the interval medians); the junction lists chr1 first although the second caller wrote chr2 first.
    expected = [["chr1", 10050, "chr2", 20050]]
    assert numpy_results["merge_callsets_translocation"] == expected
    assert fallback_results["merge_callsets_translocation"] == expected


def test_called_breakpoints_are_the_median_of_the_cluster(
    numpy_results: dict[str, Any], fallback_results: dict[str, Any]
) -> None:
    # breakpoint 1: median(1000, 1010, 1020, 1031) = 1015; breakpoint 2: median(4990, 4995, 5000, 5001) = 4997.5 -> 4997
    expected = [["chr1", 1015, 4997, "DEL", 3982]]
    assert numpy_results["sv_calls"] == expected
    assert fallback_results["sv_calls"] == expected


def test_insert_size_statistics_match_numpy(numpy_results: dict[str, Any], fallback_results: dict[str, Any]) -> None:
    _assert_rows_close(numpy_results["insert_size_stats"], _EXPECTED_INSERT_SIZE_STATS)
    _assert_rows_close(fallback_results["insert_size_stats"], _EXPECTED_INSERT_SIZE_STATS)


def test_single_insert_size_keeps_the_default_spread(fallback_results: dict[str, Any]) -> None:
    # A zero spread would collapse the discordance window to a single insert size.
    mean, std, median, mad = fallback_results["insert_size_stats"][0]
    assert (mean, std, median, mad) == (400.0, 100.0, 400.0, 0.0)


@pytest.mark.parametrize("key", ["log2_ratio", "log2_ratio_gc"])
def test_log2_ratio_matches_numpy_and_is_median_centred(
    key: str, numpy_results: dict[str, Any], fallback_results: dict[str, Any]
) -> None:
    assert fallback_results[key] == pytest.approx(numpy_results[key], rel=1e-9, abs=1e-9)
    assert statistics.median(fallback_results[key]) == pytest.approx(0.0, abs=1e-9)


def test_gc_correction_changes_the_fallback_log2_ratio(fallback_results: dict[str, Any]) -> None:
    corrected = fallback_results["log2_ratio_gc"]
    raw = fallback_results["log2_ratio"]
    assert max(abs(c - r) for c, r in zip(corrected, raw)) > 0.05


def test_coverage_track_labels_the_true_median_depth(
    numpy_results: dict[str, Any], fallback_results: dict[str, Any]
) -> None:
    # positive depths 10 and 20 -> median 15 (the upper median would print 20x)
    assert numpy_results["coverage_median_label"] == "Median: 15x"
    assert fallback_results["coverage_median_label"] == "Median: 15x"


if __name__ == "__main__":
    print(_RESULT_MARKER + json.dumps(collect_results()))
