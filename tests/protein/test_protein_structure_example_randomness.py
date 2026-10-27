"""Regression tests: the protein-structure modules do not draw from NumPy's or Python's global random state.

Audit result: ``protein/structure`` has no runtime randomness. Every use of global random state (13 in
``general.py``, 3 in ``contacts.py``, 2 in ``analysis.py``) was a docstring example such as
``coords = np.random.rand(100, 3)``. Running those examples (they are valid doctests) advanced the process-wide
NumPy generator, so a doctest run changed the draws of unrelated code, and the examples were not reproducible. The
examples now draw from a local seeded generator (``np.random.default_rng(0)``).

The checks are written for the class of problem rather than the instances: no executed example may touch the global
state, and no ``np.random.<function>`` / ``random.<function>`` call (anything but constructing a generator) may
appear anywhere in the subpackage's source.
"""

from __future__ import annotations

import doctest
import importlib
import random
import re
from pathlib import Path

import numpy as np
import pytest

import metainformant.protein.structure as structure_package

MODULES = ["analysis", "contacts", "general"]
STRUCTURE_DIR = Path(structure_package.__file__).resolve().parent

# A call on the global generators: np.random.rand(...), numpy.random.randint(...), random.sample(...), ...
# Constructing a generator (default_rng, Generator, RandomState, SeedSequence, Random) is fine.
GLOBAL_RANDOM_CALL = re.compile(
    r"\b(?:np|numpy)\.random\.(?!default_rng\b|Generator\b|SeedSequence\b|PCG64\b|RandomState\b)\w+\s*\("
    r"|(?<![\w.])random\.(?!Random\b|SystemRandom\b)\w+\s*\("
)


@pytest.mark.parametrize("module_name", MODULES)
def test_running_the_examples_leaves_the_global_random_state_alone(module_name: str) -> None:
    module = importlib.import_module(f"metainformant.protein.structure.{module_name}")
    numpy_state = np.random.get_state()
    python_state = random.getstate()

    result = doctest.testmod(module, optionflags=doctest.ELLIPSIS, verbose=False)

    assert result.attempted > 0
    assert result.failed == 0
    after = np.random.get_state()
    assert after[0] == numpy_state[0]
    assert (after[1] == numpy_state[1]).all() and after[2:] == numpy_state[2:]
    assert random.getstate() == python_state


def test_no_module_of_the_subpackage_calls_the_global_random_functions() -> None:
    offenders = []
    for path in sorted(STRUCTURE_DIR.glob("*.py")):
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if GLOBAL_RANDOM_CALL.search(line):
                offenders.append(f"{path.name}:{number}: {line.strip()}")

    assert offenders == []


def test_the_guard_pattern_recognises_global_calls_and_accepts_seeded_generators() -> None:
    """Check the guard itself: it must flag the old spellings and let the new ones through."""
    for flagged in (
        "coords = np.random.rand(100, 3)",
        ">>> mask = numpy.random.randint(0, 2, (3, 3))",
        "np.random.seed(0)",
        "x = random.sample(range(10), 3)",
        "random.shuffle(items)",
    ):
        assert GLOBAL_RANDOM_CALL.search(flagged), flagged
    for accepted in (
        "rng = np.random.default_rng(0)",
        "coords = rng.random((100, 3))",
        "gen = random.Random(7)",
        "state = np.random.RandomState(3)",
        "from numpy.random import default_rng",
    ):
        assert not GLOBAL_RANDOM_CALL.search(accepted), accepted
