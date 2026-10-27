"""``active_module_detection`` reports real per-omic contributions, not a placeholder.

Regression tests for ``metainformant.multiomics.pathways.enrichment.active_module_detection``. Every module carried
``"omic_contributions": {}`` -- a documented "placeholder dict" -- because the function only receives one combined
score per node and so has no per-omic information to report. The optional ``omic_scores`` argument now supplies it:
the contribution of an omic is the module score recomputed from that omic's p-values alone.

The oracle is the closed form ``sum(log(alpha / p))``, which equals the implementation's ``-log(p) - (-log(alpha))``
summed over the module but is computed differently.

Real computation; no test doubles.
"""

from __future__ import annotations

import math

import pytest

from metainformant.multiomics.pathways.enrichment import active_module_detection

ALPHA = 0.05
NETWORK = {"G1": ["G2", "G3"], "G2": ["G1"], "G3": ["G1"], "G4": ["G5"], "G5": ["G4"]}
COMBINED = {"G1": 0.001, "G2": 0.002, "G3": 0.003, "G4": 0.9, "G5": 0.8}


def _module(omic_scores: dict[str, dict[str, float]] | None = None) -> dict:
    """The single module of NETWORK (G1-G3), with per-omic data only when it is given."""
    if omic_scores is None:
        modules = active_module_detection(NETWORK, COMBINED, alpha=ALPHA, n_permutations=50)
    else:
        modules = active_module_detection(NETWORK, COMBINED, alpha=ALPHA, n_permutations=50, omic_scores=omic_scores)
    (module,) = modules
    return module


class TestOmicContributions:
    def test_each_omic_contributes_its_own_module_score(self) -> None:
        omic_scores = {
            "rna": {"G1": 0.001, "G2": 0.002, "G3": 0.003},
            "protein": {"G1": 0.04, "G2": 0.5},  # G3 was not measured: it counts as p = 1, as in `scores`
        }

        contributions = _module(omic_scores=omic_scores)["omic_contributions"]

        assert contributions == pytest.approx(
            {
                "rna": math.log(ALPHA / 0.001) + math.log(ALPHA / 0.002) + math.log(ALPHA / 0.003),
                "protein": math.log(ALPHA / 0.04) + math.log(ALPHA / 0.5) + math.log(ALPHA / 1.0),
            },
            rel=1e-12,
        )

    def test_the_driving_omic_scores_above_zero_and_a_non_supporting_omic_below(self) -> None:
        omic_scores = {
            "rna": {"G1": 1e-6, "G2": 1e-5, "G3": 1e-4},
            "methylation": {"G1": 0.7, "G2": 0.8, "G3": 0.9},
        }

        contributions = _module(omic_scores=omic_scores)["omic_contributions"]

        assert contributions["rna"] > 0
        assert contributions["methylation"] < 0
        assert contributions["rna"] > contributions["methylation"]

    def test_a_single_omic_equal_to_the_combined_scores_reproduces_the_module_score(self) -> None:
        module = _module(omic_scores={"all": dict(COMBINED)})

        assert module["omic_contributions"]["all"] == pytest.approx(module["module_score"], rel=1e-12)

    def test_the_contributions_cover_exactly_the_omics_that_were_given(self) -> None:
        omic_scores = {"a": {"G1": 0.01}, "b": {}, "c": {"G2": 0.01, "G3": 0.01}}

        assert set(_module(omic_scores=omic_scores)["omic_contributions"]) == {"a", "b", "c"}

    def test_per_omic_scores_do_not_change_which_modules_are_found_or_how_they_score(self) -> None:
        plain = active_module_detection(NETWORK, COMBINED, alpha=ALPHA, n_permutations=50)
        with_omics = active_module_detection(
            NETWORK, COMBINED, alpha=ALPHA, n_permutations=50, omic_scores={"rna": {"G1": 0.001}}
        )

        assert [(m["module_genes"], m["module_score"], m["p_value"]) for m in with_omics] == [
            (m["module_genes"], m["module_score"], m["p_value"]) for m in plain
        ]

    def test_without_per_omic_scores_the_field_is_empty_because_nothing_can_be_attributed(self) -> None:
        assert _module()["omic_contributions"] == {}
        assert _module(omic_scores={})["omic_contributions"] == {}

    def test_every_module_gets_its_own_contributions(self) -> None:
        network = {"A": ["B"], "B": ["A"], "X": ["Y"], "Y": ["X"]}
        combined = {"A": 0.001, "B": 0.001, "X": 0.002, "Y": 0.002}
        omic_scores = {
            "rna": {"A": 0.001, "B": 0.001, "X": 0.5, "Y": 0.5},
            "prot": {"A": 0.5, "B": 0.5, "X": 0.002, "Y": 0.002},
        }

        modules = active_module_detection(network, combined, alpha=ALPHA, n_permutations=50, omic_scores=omic_scores)

        by_genes = {tuple(m["module_genes"]): m["omic_contributions"] for m in modules}
        assert by_genes[("A", "B")]["rna"] == pytest.approx(2 * math.log(ALPHA / 0.001))
        assert by_genes[("A", "B")]["prot"] == pytest.approx(2 * math.log(ALPHA / 0.5))
        assert by_genes[("X", "Y")]["rna"] == pytest.approx(2 * math.log(ALPHA / 0.5))
        assert by_genes[("X", "Y")]["prot"] == pytest.approx(2 * math.log(ALPHA / 0.002))
