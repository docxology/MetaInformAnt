"""Regression tests for the protein sequence-logo letter heights.

A logo column is ``log2(20) - H`` bits tall, where ``H`` is the Shannon entropy of
the residue frequencies, and each residue's letter is its frequency times that
height.  The heights used to be ``log2(20) + H``, so the *more variable* a column
was the taller it was drawn.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from metainformant.protein.visualization.general import plot_sequence_logo  # noqa: E402

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
MAX_BITS = math.log2(len(AMINO_ACIDS))


def _letter_heights(ax) -> dict[str, np.ndarray]:
    """Bar heights per amino acid (one BarContainer per residue, one bar per position)."""
    assert len(ax.containers) == len(AMINO_ACIDS)
    return {
        aa: np.array([bar.get_height() for bar in container.patches])
        for aa, container in zip(AMINO_ACIDS, ax.containers)
    }


def _column_heights(ax) -> np.ndarray:
    return np.sum(list(_letter_heights(ax).values()), axis=0)


def test_fully_conserved_column_is_log2_20_bits_tall():
    ax = plot_sequence_logo(["AAAA"] * 4)

    assert _column_heights(ax) == pytest.approx([MAX_BITS] * 4)
    assert _letter_heights(ax)["A"] == pytest.approx([MAX_BITS] * 4)


def test_two_equally_frequent_residues_lose_exactly_one_bit():
    ax = plot_sequence_logo(["A", "A", "C", "C"])

    heights = _letter_heights(ax)
    assert _column_heights(ax) == pytest.approx([MAX_BITS - 1.0])
    assert heights["A"] == pytest.approx([0.5 * (MAX_BITS - 1.0)])
    assert heights["C"] == pytest.approx([0.5 * (MAX_BITS - 1.0)])


def test_column_using_all_twenty_residues_equally_carries_no_information():
    ax = plot_sequence_logo(list(AMINO_ACIDS))

    assert _column_heights(ax) == pytest.approx([0.0], abs=1e-9)


def test_variable_columns_are_shorter_than_conserved_columns():
    # column 0: A x4 (H = 0), column 1: A/C/D/E once each (H = 2 bits), column 2: A x3 + C (H = 0.811 bits)
    ax = plot_sequence_logo(["AAA", "ACA", "ADA", "AEC"])

    heights = _column_heights(ax)
    entropy_column_2 = -(0.75 * math.log2(0.75) + 0.25 * math.log2(0.25))
    assert heights == pytest.approx([MAX_BITS, MAX_BITS - 2.0, MAX_BITS - entropy_column_2])
    assert heights[0] > heights[2] > heights[1]


def test_total_height_never_exceeds_the_maximum_and_never_goes_negative():
    rng = np.random.default_rng(3)
    sequences = ["".join(rng.choice(list(AMINO_ACIDS), size=12)) for _ in range(40)]

    heights = _column_heights(plot_sequence_logo(sequences))

    assert np.all(heights >= -1e-9)
    assert np.all(heights <= MAX_BITS + 1e-9)


def test_logo_is_saved_when_an_output_path_is_given(tmp_path):
    out = tmp_path / "logo.png"

    plot_sequence_logo(["ACDE", "ACDF"], output_path=out)

    assert out.exists() and out.stat().st_size > 0
