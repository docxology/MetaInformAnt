"""A ``tf_list`` restricts the regulators of the inferred GRN; it never silently widens to every gene.

``infer_grn_correlation``, ``infer_grn_mutual_info`` and ``infer_grn_regression`` replaced a ``tf_list`` that matched
no gene in ``gene_names`` (a misspelt or wrongly-cased catalogue, or an empty list) with *all* genes. The caller asked
for a regulator-restricted network and got the full all-against-all one, the correlation function with a log line and
the other two with nothing at all. Repeated entries of ``tf_list`` also produced duplicate edges (and, for regression,
duplicate coefficient columns).

A list that matches nothing now raises ``ValueError``; names that are not in ``gene_names`` are still ignored when at
least one regulator matches (a genome-wide TF catalogue is longer than the measured genes) and repeated names count
once. Real computation on planted expression data; no test doubles.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.networks.regulatory.grn_inference import (
    infer_grn_correlation,
    infer_grn_mutual_info,
    infer_grn_regression,
)

NAMES = ["G0", "G1", "G2", "G3"]


def _expression() -> np.ndarray:
    rng = np.random.default_rng(0)
    base = rng.normal(size=40)
    return np.vstack([base, base + 0.1 * rng.normal(size=40), -base + 0.1 * rng.normal(size=40), rng.normal(size=40)])


INFERENCE = {
    "correlation": lambda tf_list: infer_grn_correlation(_expression(), NAMES, tf_list=tf_list, threshold=0.5),
    "mutual_info": lambda tf_list: infer_grn_mutual_info(_expression(), NAMES, tf_list=tf_list),
    "regression": lambda tf_list: infer_grn_regression(_expression(), NAMES, tf_list=tf_list, alpha=0.01),
}


@pytest.mark.parametrize("method", sorted(INFERENCE))
@pytest.mark.parametrize("tf_list", [["NOT_A_GENE"], ["g0", "g1"], [], ["x", "y", "z", "w", "v", "u", "t"]])
def test_a_tf_list_that_matches_no_gene_raises(method: str, tf_list: list[str]) -> None:
    with pytest.raises(ValueError, match="tf_list"):
        INFERENCE[method](tf_list)


def test_the_error_names_what_was_asked_for() -> None:
    with pytest.raises(ValueError, match="'NOT_A_GENE'"):
        INFERENCE["correlation"](["NOT_A_GENE"])
    with pytest.raises(ValueError, match="empty list"):
        INFERENCE["mutual_info"]([])


@pytest.mark.parametrize("method", sorted(INFERENCE))
def test_unknown_names_are_ignored_when_one_regulator_matches(method: str) -> None:
    result = INFERENCE[method](["G0", "NOT_A_GENE"])

    assert result["edges"], "G0 is strongly related to G1 and G2, so it must have edges"
    assert {edge["source"] for edge in result["edges"]} == {"G0"}


@pytest.mark.parametrize("method", sorted(INFERENCE))
def test_repeated_names_do_not_duplicate_edges(method: str) -> None:
    once = INFERENCE[method](["G0"])
    repeated = INFERENCE[method](["G0", "G0", "G0"])

    assert repeated["edges"] == once["edges"]
    assert repeated["n_edges"] == once["n_edges"]


@pytest.mark.parametrize("method", sorted(INFERENCE))
def test_no_tf_list_still_means_every_gene_is_a_regulator(method: str) -> None:
    everything = {
        "correlation": infer_grn_correlation(_expression(), NAMES, threshold=0.5),
        "mutual_info": infer_grn_mutual_info(_expression(), NAMES),
        "regression": infer_grn_regression(_expression(), NAMES, alpha=0.01),
    }[method]

    assert len({edge["source"] for edge in everything["edges"]}) > 1


def test_regression_coefficient_columns_follow_the_unique_regulators() -> None:
    result = INFERENCE["regression"](["G0", "G0", "G1"])

    # one coefficient per unique regulator (G0, G1) for each of the four target genes
    assert all(len(row) == 2 for row in result["coefficients"])
    assert len(result["coefficients"]) == len(NAMES)
