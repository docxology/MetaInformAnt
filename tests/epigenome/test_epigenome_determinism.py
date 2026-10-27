"""Reproducibility of the epigenome analyses: no unseeded randomness, no hash-order dependence.

Audit of ``src/metainformant/epigenome`` for random number use: nothing draws from the global ``random`` module or from
``numpy.random``. The one generator (the k-means++ start of ``learn_chromatin_states``) is a private
``random.Random(42)``, so it is seeded and does not touch global state. ``TestNoUnseededRandomness`` keeps it that way.

The audit did find results whose *order* depended on the iteration order of a ``set`` of chromosome names, which changes
from one Python process to the next (string hashing is randomized per process): ``find_overlapping_peaks`` (order of the
returned pairs), ``find_differentially_methylated_regions`` (order of the regions) and ``merge_tracks`` (order of the
chromosomes of the merged track). ``TestResultsDoNotDependOnTheHashSeed`` runs those three in fresh interpreters with
different ``PYTHONHASHSEED`` values and requires identical, sorted output.

REAL IMPLEMENTATION.
"""

from __future__ import annotations

import ast
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List

import pytest

import metainformant
import metainformant.epigenome

# ---------------------------------------------------------------------------
# Static audit
# ---------------------------------------------------------------------------


def _epigenome_sources() -> List[Path]:
    return sorted(Path(metainformant.epigenome.__file__).parent.rglob("*.py"))


def _rng_findings(path: Path) -> List[str]:
    """Uses of global or unseeded random number generation in one source file."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    findings: List[str] = []

    for node in ast.walk(tree):
        # import numpy.random / from numpy import random / from numpy.random import x
        if isinstance(node, ast.Import):
            findings.extend(
                f"{path.name}:{node.lineno} import {a.name}" for a in node.names if a.name.startswith("numpy.random")
            )
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if module.startswith("numpy.random") or (module == "numpy" and any(a.name == "random" for a in node.names)):
                findings.append(f"{path.name}:{node.lineno} from {module} import ...")
            if module == "random":
                # only the Random class (an explicit generator) may be imported; the functions use the global generator
                findings.extend(
                    f"{path.name}:{node.lineno} from random import {a.name}" for a in node.names if a.name != "Random"
                )

        elif isinstance(node, ast.Attribute):
            owner = node.value
            # random.shuffle, random.random, random.seed, ...: the global generator
            if isinstance(owner, ast.Name) and owner.id == "random" and node.attr != "Random":
                findings.append(f"{path.name}:{node.lineno} random.{node.attr}")
            # np.random.* / numpy.random.*
            if (
                isinstance(owner, ast.Attribute)
                and owner.attr == "random"
                and isinstance(owner.value, ast.Name)
                and owner.value.id in {"np", "numpy"}
            ):
                findings.append(f"{path.name}:{node.lineno} {owner.value.id}.random.{node.attr}")

        elif isinstance(node, ast.Call):
            func = node.func
            # random.Random() with no seed, numpy default_rng() with no seed
            if (
                isinstance(func, ast.Attribute)
                and func.attr in {"Random", "default_rng", "RandomState"}
                and not node.args
                and not node.keywords
            ):
                findings.append(f"{path.name}:{node.lineno} unseeded {func.attr}()")

    return findings


class TestNoUnseededRandomness:
    def test_the_package_has_python_sources_to_audit(self) -> None:
        names = {p.name for p in _epigenome_sources()}

        assert {"chipseq.py", "atacseq.py", "methylation.py", "peak_detection.py", "tracks.py", "workflow.py"} <= names

    def test_no_global_or_unseeded_random_number_generation(self) -> None:
        findings = [finding for path in _epigenome_sources() for finding in _rng_findings(path)]

        assert findings == []

    def test_the_audit_recognises_what_it_forbids(self, tmp_path: Path) -> None:
        # A file that does all the forbidden things must be flagged, or the audit above would pass vacuously
        offender = tmp_path / "offender.py"
        offender.write_text(
            "import random\n"
            "import numpy as np\n"
            "from random import shuffle\n"
            "from numpy.random import default_rng\n"
            "x = random.random()\n"
            "y = np.random.rand(3)\n"
            "z = random.Random()\n"
            "ok = random.Random(42)\n"
        )

        findings = _rng_findings(offender)

        # lines 3 to 7; the seeded generator on line 8 is fine
        assert sorted(int(finding.split(":")[1].split(" ")[0]) for finding in findings) == [3, 4, 5, 6, 7]


# ---------------------------------------------------------------------------
# Hash-seed independence
# ---------------------------------------------------------------------------

# Six chromosome names whose set order varies with the hash seed. They are deliberately not in sorted order.
_CHROMOSOMES = ["chrX", "chr2", "chr10", "chr1", "chrY", "chr3"]

_SCENARIO = f"""
import json, logging
logging.disable(logging.CRITICAL)
from metainformant.epigenome.assays.chipseq import ChIPPeak, find_overlapping_peaks
from metainformant.epigenome.assays.methylation import MethylationSite, find_differentially_methylated_regions
from metainformant.epigenome.analysis.tracks import GenomicTrack, merge_tracks

chromosomes = {_CHROMOSOMES!r}

peaks1 = [ChIPPeak(c, 100, 200) for c in chromosomes]
peaks2 = [ChIPPeak(c, 150, 250) for c in chromosomes]
pairs = [[a.chromosome, a.start, b.start] for a, b in find_overlapping_peaks(peaks1, peaks2)]

def sites(methylated):
    return {{c: [MethylationSite(c, 1000 + 10 * i, methylated, 10) for i in range(5)] for c in chromosomes}}
dmrs = find_differentially_methylated_regions(sites(0), sites(10), min_sites=3)
regions = [[d["chromosome"], d["start"]] for d in dmrs]

def track():
    t = GenomicTrack()
    for c in chromosomes:
        t.add_feature(c, 0, 100, value=1.0)
    return t
union = merge_tracks([track(), track()], operation="union").get_chromosomes()
intersection = merge_tracks([track(), track()], operation="intersection").get_chromosomes()

print(json.dumps({{"pairs": pairs, "regions": regions, "union": union, "intersection": intersection}}))
"""


def _run_scenario(hash_seed: int) -> Dict[str, Any]:
    """Run the scenario in a fresh interpreter that imports the same metainformant as this test process."""
    package_root = Path(metainformant.__file__).resolve().parents[1]
    env = {
        "PATH": os.environ.get("PATH", ""),
        "HOME": os.environ.get("HOME", ""),
        "TMPDIR": os.environ.get("TMPDIR", tempfile.gettempdir()),
        "PYTHONHASHSEED": str(hash_seed),
        "PYTHONPATH": str(package_root),
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    completed = subprocess.run(
        [sys.executable, "-c", _SCENARIO], capture_output=True, text=True, check=True, env=env, timeout=300
    )
    parsed: Dict[str, Any] = json.loads(completed.stdout.strip().splitlines()[-1])
    return parsed


class TestResultsDoNotDependOnTheHashSeed:
    @pytest.mark.slow
    def test_output_is_identical_and_ordered_whatever_the_hash_seed(self) -> None:
        outputs = [_run_scenario(seed) for seed in (1, 2, 3)]

        assert all(output == outputs[0] for output in outputs[1:])

        output = outputs[0]
        # overlapping pairs: chromosomes in the order they first appear in the first peak set
        assert [pair[0] for pair in output["pairs"]] == _CHROMOSOMES
        # DMRs and merged-track chromosomes: sorted by name
        assert [region[0] for region in output["regions"]] == sorted(_CHROMOSOMES)
        assert output["union"] == sorted(_CHROMOSOMES)
        assert output["intersection"] == sorted(_CHROMOSOMES)
