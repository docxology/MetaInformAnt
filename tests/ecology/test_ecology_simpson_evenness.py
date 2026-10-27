"""Simpson's evenness is (1 / sum p_i^2) / S, so a perfectly even community scores 1."""

from __future__ import annotations

import pytest

from metainformant.ecology.analysis.community import calculate_evenness


def test_even_community_scores_one() -> None:
    assert calculate_evenness([1, 1, 1, 1], "simpson") == pytest.approx(1.0)


def test_uneven_community_matches_hand_computation() -> None:
    # p = 0.8, 0.1, 0.1: sum p^2 = 0.66, so E = (1 / 0.66) / 3.
    assert calculate_evenness([8, 1, 1], "simpson") == pytest.approx(1 / 0.66 / 3)
