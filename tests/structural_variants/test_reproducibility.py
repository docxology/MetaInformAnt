"""Reproducibility audit of the structural-variant analysis paths.

Two kinds of run-to-run variation were looked for:

* Random number generators. The package draws no random numbers at all (the CBS p-value is analytic, the
  merging and clustering are greedy and ordered), and ``test_no_random_number_generator_in_the_package`` keeps
  it that way: an analysis path that needs randomness must take an explicit seed, and this list is then
  updated on purpose.
* Hash order. Python randomises the order of a ``set`` of strings per process (``PYTHONHASHSEED``). Four
  results were built with ``list(set(...))`` and so came out in a different order on every run:
  ``MergedVariant.caller_names`` (``merge_callsets`` and ``survivor_merge``), ``SVEvidence.evidence_reads`` and
  ``TADPrediction.genes_in_affected_tads``. They are sorted now.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

import metainformant.structural_variants as structural_variants

# Everything below runs in a fresh interpreter, so the hash seed can be chosen.
SCRIPT = r"""
import json
from dataclasses import asdict

from metainformant.structural_variants.annotation.functional_impact import predict_tad_disruption
from metainformant.structural_variants.detection.sv_calling import call_structural_variants
from metainformant.structural_variants.filtering.merge import merge_callsets, survivor_merge

callers = ["delly", "manta", "lumpy", "gridss", "svaba", "pindel"]
callsets = {
    name: [{"chrom": "chr1", "start": 1000 + 10 * k, "end": 5000 + 10 * k, "sv_type": "DEL", "_caller": name}]
    for k, name in enumerate(callers)
}
merged, _ = merge_callsets(callsets)
survivors, _ = survivor_merge([v[0] for v in callsets.values()], max_distance=500, min_callers=2)

# One translocation per caller, half of them written from the other side.
def junction(name, k):
    if k % 2 == 0:
        first, second = ("chr1", 1000 + 10 * k), ("chr2", 5000 + 10 * k)
    else:
        first, second = ("chr2", 5000 + 10 * k), ("chr1", 1000 + 10 * k)
    return {
        "chrom": first[0], "start": first[1], "end": first[1], "sv_type": "TRA",
        "chrom2": second[0], "pos2": second[1], "_caller": name,
    }


junctions = {name: [junction(name, k)] for k, name in enumerate(callers)}
junction_merged, _ = merge_callsets(junctions)
junction_survivors, _ = survivor_merge([v[0] for v in junctions.values()], max_distance=500, min_callers=2)

reads = [
    {
        "name": f"read_{i:02d}", "chrom": "chr1", "pos": 1000 + i, "mapq": 60, "cigar": "50M60S",
        "sa_tag": "chr1,5000,-,60S50M,60,0;", "insert_size": 400, "mate_chrom": "chr1", "mate_pos": 1400,
        "is_reverse": False, "mate_is_reverse": True,
    }
    for i in range(12)
]
calls = call_structural_variants(reads, min_support=3)

tad = predict_tad_disruption(
    {"chrom": "chr1", "start": 50_000, "end": 150_000, "sv_type": "DEL"},
    [{"chrom": "chr1", "start": 99_000, "end": 101_000, "genes": ["BRCA1", "TP53", "EGFR", "MYC", "KRAS", "PTEN"]}],
)

print(json.dumps({
    "merge_callsets": [m.caller_names for m in merged],
    "survivor_merge": [m.caller_names for m in survivors],
    "junction_merge": [[m.chrom, m.start, m.chrom2, m.pos2, m.caller_names] for m in junction_merged],
    "junction_survivor_merge": [[m.chrom, m.start, m.chrom2, m.pos2, m.caller_names] for m in junction_survivors],
    "evidence_reads": [c.evidence.evidence_reads for c in calls],
    "calls": [json.dumps(asdict(c), sort_keys=True, default=str) for c in calls],
    "tad_genes": tad.genes_in_affected_tads,
}))
"""


def _start_with_hash_seed(seed: int, workdir: Path) -> subprocess.Popen[str]:
    # The child must import the very copy of the package that this process is testing: the directory that
    # holds ``metainformant`` (``src`` in the repository), not whichever copy the environment lists first.
    package_root = Path(structural_variants.__file__).parents[2]
    env = dict(os.environ, PYTHONHASHSEED=str(seed), PYTHONPATH=str(package_root))
    return subprocess.Popen(
        [sys.executable, "-c", SCRIPT], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=workdir, env=env
    )


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> list[dict[str, list]]:
    workdir = tmp_path_factory.mktemp("hash_seeds")
    children = [_start_with_hash_seed(seed, workdir) for seed in (0, 1, 2, 3)]  # in parallel
    results: list[dict[str, list]] = []
    for child in children:
        stdout, stderr = child.communicate(timeout=300)
        assert child.returncode == 0, stderr[-3000:]
        results.append(json.loads(stdout.strip().splitlines()[-1]))
    return results


class TestHashOrder:
    def test_results_are_identical_whatever_the_hash_seed(self, runs: list[dict[str, list]]) -> None:
        for later in runs[1:]:
            assert later == runs[0]

    def test_unordered_collections_come_out_sorted(self, runs: list[dict[str, list]]) -> None:
        first = runs[0]
        assert first["merge_callsets"] == [["delly", "gridss", "lumpy", "manta", "pindel", "svaba"]]
        assert first["survivor_merge"] == [["delly", "gridss", "lumpy", "manta", "pindel", "svaba"]]
        # Six callers, first breakpoints 1000..1050 and second breakpoints 5000..5050 whichever side each wrote
        # first: the medians are 1025 and 5025.
        everyone = ["delly", "gridss", "lumpy", "manta", "pindel", "svaba"]
        assert first["junction_merge"] == [["chr1", 1025, "chr2", 5025, everyone]]
        assert first["junction_survivor_merge"] == [["chr1", 1025, "chr2", 5025, everyone]]
        assert first["evidence_reads"] == [[f"read_{i:02d}" for i in range(12)]]
        assert first["tad_genes"] == ["BRCA1", "EGFR", "KRAS", "MYC", "PTEN", "TP53"]

    def test_the_comparison_is_not_trivial(self, runs: list[dict[str, list]]) -> None:
        assert runs[0]["calls"]  # a call was made, so evidence_reads was compared


# ---------------------------------------------------------------------------
# Random number generators
# ---------------------------------------------------------------------------

# Anything that draws random numbers, or an ambient source of variation.
RANDOMNESS = re.compile(
    r"np\.random|numpy\.random|default_rng|RandomState|^\s*import random\b|^\s*from random\b"
    r"|\bsecrets\b|\buuid\b|os\.urandom|time\.time\(|datetime\.now|\bhash\(",
    re.MULTILINE,
)


def test_no_random_number_generator_in_the_package() -> None:
    package = Path(structural_variants.__file__).parent
    sources = sorted(path for path in package.rglob("*.py") if "__pycache__" not in path.parts)
    assert len(sources) >= 14  # every module of the package is scanned

    offenders = {
        str(path.relative_to(package)): [m.group(0).strip() for m in RANDOMNESS.finditer(path.read_text())]
        for path in sources
        if RANDOMNESS.search(path.read_text())
    }
    assert offenders == {}
