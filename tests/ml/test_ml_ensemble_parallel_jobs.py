"""Regression tests: the ensemble in ``cross_validate_biological`` follows ``parallel_jobs``.

``cross_validate_biological(method="ensemble")`` trained its voting ensemble with ``n_jobs=-1`` unless an ``n_jobs``
keyword was passed, whatever ``parallel_jobs`` said. ``parallel_jobs=1`` ("use 1 in restricted or resource-sensitive
environments") therefore still started a joblib process pool for the ensemble, which fails outright where process
pools are unavailable (such as a sandbox that forbids ``os.sysconf("SC_SEM_NSEMS_MAX")``) and oversubscribes the
machine elsewhere. The ensemble now uses ``parallel_jobs`` unless ``n_jobs`` is given explicitly.

Where process pools are unavailable the in-process check fails on the old code; elsewhere the subprocess check does,
because it asserts that no worker process was started. The oracle is a hand-assembled soft-voting ensemble with the
same members run sequentially through scikit-learn's own ``cross_val_score``. No test doubles.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from sklearn.datasets import make_classification
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier, VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score

import metainformant
from metainformant.ml.models.classification import cross_validate_biological

X, Y = make_classification(n_samples=90, n_features=5, n_informative=3, n_redundant=0, random_state=0)
SEED = 3


def _oracle_scores() -> list[float]:
    ensemble = VotingClassifier(
        estimators=[
            ("rf", RandomForestClassifier(n_estimators=10, random_state=SEED)),
            ("gb", GradientBoostingClassifier(n_estimators=10, random_state=SEED)),
            ("lr", LogisticRegression(random_state=SEED, max_iter=1000)),
        ],
        voting="soft",
        n_jobs=1,
    )
    folds = StratifiedKFold(n_splits=3, shuffle=True, random_state=SEED)
    return cross_val_score(ensemble, X, Y, cv=folds, scoring="accuracy", n_jobs=1).tolist()


@pytest.mark.parametrize("parallel_jobs", [1, None])
def test_the_ensemble_runs_sequentially_when_parallel_jobs_asks_for_it(parallel_jobs: int | None) -> None:
    result = cross_validate_biological(
        X, Y, method="ensemble", cv_folds=3, random_state=SEED, parallel_jobs=parallel_jobs
    )  # no n_jobs keyword: the ensemble has to follow parallel_jobs on its own

    assert result["cross_validation"]["accuracy"]["scores"] == pytest.approx(_oracle_scores())


SCRIPT = """
import json, multiprocessing

from sklearn.datasets import make_classification

from metainformant.ml.models.classification import cross_validate_biological

X, y = make_classification(n_samples=90, n_features=5, n_informative=3, n_redundant=0, random_state=0)
result = cross_validate_biological(X, y, method="ensemble", cv_folds=3, random_state=3, parallel_jobs=1)
print(json.dumps({"scores": result["cross_validation"]["accuracy"]["scores"],
                  "worker_processes": len(multiprocessing.active_children())}))
"""


def test_parallel_jobs_one_starts_no_worker_processes() -> None:
    package_root = str(Path(metainformant.__file__).parent.parent)
    existing = os.environ.get("PYTHONPATH")
    env = {**os.environ, "PYTHONPATH": os.pathsep.join([package_root] + ([existing] if existing else []))}

    completed = subprocess.run(
        [sys.executable, "-c", SCRIPT], capture_output=True, text=True, env=env, timeout=300, check=False
    )

    assert completed.returncode == 0, completed.stderr
    observed = json.loads(completed.stdout.strip().splitlines()[-1])
    assert observed["worker_processes"] == 0
    assert observed["scores"] == pytest.approx(_oracle_scores())


def test_an_explicit_n_jobs_keyword_still_wins_over_parallel_jobs() -> None:
    result = cross_validate_biological(
        X, Y, method="ensemble", cv_folds=3, random_state=SEED, parallel_jobs=1, n_jobs=1
    )

    assert result["cross_validation"]["accuracy"]["scores"] == pytest.approx(_oracle_scores())
