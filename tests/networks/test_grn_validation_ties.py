"""``validate_grn`` AUROC and AUPRC do not depend on the order of tied scores (or on ``PYTHONHASHSEED``).

The candidate edges of ``validate_grn`` are the ordered pairs of a ``set`` of node names, scored with the predicted edge
weight or 0.0, then sorted by score. The curve code walked the sorted items one by one, so the order of tied items
(the set's iteration order, which changes with the string hash seed) decided the area. An unweighted prediction scores
every predicted edge 1.0 and every other pair 0.0: two tie groups, and an AUROC that ranged from 0.74 to 0.96 and an
AUPRC from 0.50 to 0.83 over eight hash seeds when the answer is 0.8148 and 0.4778.

The curves now move between distinct scores only (a straight ROC segment, one precision evaluation per threshold), which
is the Mann-Whitney AUC with ties counted half and scikit-learn's step-wise average precision. The oracles are
``sklearn.metrics.roc_auc_score`` / ``average_precision_score`` on the same labels and scores. The hash-seed cases run
the real module in subprocesses with a pinned path and ``PYTHONHASHSEED``. No test doubles.
"""

from __future__ import annotations

import json
import os
import random
import subprocess
import sys
from pathlib import Path

import pytest

sklearn_metrics = pytest.importorskip("sklearn.metrics")

import metainformant  # noqa: E402
from metainformant.networks.regulatory.grn_inference import _compute_auprc, _compute_auroc, validate_grn  # noqa: E402

GENES = [f"G{i}" for i in range(6)]
KNOWN = [("G0", "G1"), ("G0", "G2"), ("G3", "G4")]
PREDICTED = [("G0", "G1"), ("G0", "G2"), ("G4", "G5")]  # unweighted: two of three are right


def _edges(pairs, weights=None):
    if weights is None:
        return [{"source": s, "target": t} for s, t in pairs]
    return [{"source": s, "target": t, "weight": w} for (s, t), w in zip(pairs, weights)]


def _oracle(predicted: dict, known: set, genes: list[str]) -> tuple[float, float]:
    pairs = [(a, b) for a in genes for b in genes if a != b]
    labels = [pair in known for pair in pairs]
    scores = [abs(predicted.get(pair, 0.0)) for pair in pairs]
    return (
        float(sklearn_metrics.roc_auc_score(labels, scores)),
        float(sklearn_metrics.average_precision_score(labels, scores)),
    )


def _subprocess_env(hash_seed: int) -> dict[str, str]:
    package_root = str(Path(metainformant.__file__).parent.parent)
    existing = os.environ.get("PYTHONPATH")
    return {
        **os.environ,
        "PYTHONHASHSEED": str(hash_seed),
        "PYTHONPATH": os.pathsep.join([package_root] + ([existing] if existing else [])),
    }


HASH_SEED_SCRIPT = """
import json, sys

from metainformant.networks.regulatory.grn_inference import validate_grn

payload = json.load(sys.stdin)
print(json.dumps(validate_grn(payload["predicted"], payload["known"])))
"""


def test_an_unweighted_prediction_has_the_same_scores_under_every_hash_seed() -> None:
    expected_auroc, expected_auprc = _oracle({pair: 1.0 for pair in PREDICTED}, set(KNOWN), GENES)
    assert expected_auroc == pytest.approx(0.814815, abs=1e-6) and expected_auprc == pytest.approx(0.477778, abs=1e-6)
    payload = json.dumps({"predicted": _edges(PREDICTED), "known": _edges(KNOWN)})

    for hash_seed in range(8):
        completed = subprocess.run(
            [sys.executable, "-c", HASH_SEED_SCRIPT],
            input=payload,
            capture_output=True,
            text=True,
            env=_subprocess_env(hash_seed),
            timeout=300,
            check=False,
        )
        assert completed.returncode == 0, completed.stderr
        result = json.loads(completed.stdout.strip().splitlines()[-1])

        assert result["auroc"] == pytest.approx(expected_auroc, abs=1e-6), hash_seed
        assert result["auprc"] == pytest.approx(expected_auprc, abs=1e-6), hash_seed
        assert result["precision"] == pytest.approx(2 / 3, abs=1e-6) and result["recall"] == pytest.approx(
            2 / 3, abs=1e-6
        )


@pytest.mark.parametrize("seed", range(6))
def test_validate_grn_matches_scikit_learn_for_weights_with_many_ties(seed: int) -> None:
    rng = random.Random(seed)
    genes = [f"g{i}" for i in range(9)]
    pairs = [(a, b) for a in genes for b in genes if a != b]
    known = set(rng.sample(pairs, 14))
    predicted_pairs = rng.sample(pairs, 20)
    weights = [rng.choice([0.2, 0.5, 0.5, 1.0, -1.0]) for _ in predicted_pairs]  # |weight| is the score
    # every node must appear in the edge lists so the candidate universe is the full set of ordered pairs
    present = {n for pair in predicted_pairs + sorted(known) for n in pair}
    assert present == set(genes)

    result = validate_grn(_edges(predicted_pairs, weights), _edges(sorted(known)))

    auroc, auprc = _oracle({pair: w for pair, w in zip(predicted_pairs, weights)}, known, genes)
    assert result["auroc"] == pytest.approx(auroc, abs=1e-6)
    assert result["auprc"] == pytest.approx(auprc, abs=1e-6)


@pytest.mark.parametrize("seed", range(8))
def test_the_curve_areas_do_not_depend_on_the_order_of_tied_items(seed: int) -> None:
    rng = random.Random(100 + seed)
    scored = [(float(rng.randint(0, 3)), rng.random() < 0.35) for _ in range(rng.randint(30, 120))]
    labels = [positive for _, positive in scored]
    scores = [score for score, _ in scored]
    if not any(labels) or all(labels):
        pytest.skip("degenerate draw")

    shuffled_a = sorted(scored, key=lambda item: item[0], reverse=True)
    shuffled_b = list(shuffled_a)
    for items in (shuffled_a, shuffled_b):  # reorder inside every tie group, differently
        rng.shuffle(items)
        items.sort(key=lambda item: item[0], reverse=True)

    expected_auroc = float(sklearn_metrics.roc_auc_score(labels, scores))
    expected_auprc = float(sklearn_metrics.average_precision_score(labels, scores))
    for ordered in (shuffled_a, shuffled_b):
        assert _compute_auroc(ordered) == pytest.approx(expected_auroc, abs=1e-12)
        assert _compute_auprc(ordered) == pytest.approx(expected_auprc, abs=1e-12)


def test_one_class_inputs_keep_their_conventions() -> None:
    all_positive = [(1.0, True), (0.5, True)]
    all_negative = [(1.0, False), (0.5, False)]

    assert _compute_auroc(all_positive) == 0.5 and _compute_auroc(all_negative) == 0.5
    assert _compute_auprc(all_negative) == 0.0
    assert _compute_auprc(all_positive) == pytest.approx(1.0)


def test_a_perfect_prediction_still_scores_one() -> None:
    result = validate_grn(_edges(KNOWN), _edges(KNOWN))

    assert result["auroc"] == pytest.approx(1.0) and result["auprc"] == pytest.approx(1.0)
    assert result["precision"] == 1.0 and result["recall"] == 1.0 and result["f1"] == 1.0
