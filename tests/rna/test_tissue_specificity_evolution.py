"""Tests for descriptive tissue-specificity evolution (Mantica 2024 adaptation)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis import tissue_specificity_evolution as tse


def _synthetic_expression() -> pd.DataFrame:
    """Tiny deterministic expression matrix (orthogroups x tissues)."""
    return pd.DataFrame(
        {
            "brain": [10.0, 0.0, 0.0, 0.0],
            "muscle": [0.0, 10.0, 0.0, 1.0],
            "gut": [0.0, 0.0, 10.0, 1.0],
            "ovary": [0.0, 0.0, 1.0, 1.0],
        },
        index=["og_brain", "og_muscle", "og_gut", "og_housekeeping"],
    )


class TestComputeTau:
    def test_perfectly_specific_gene_has_tau_one(self) -> None:
        tau = tse.compute_tau(_synthetic_expression())
        assert tau.loc["og_brain", "tau"] == pytest.approx(1.0)
        assert tau.loc["og_muscle", "tau"] == pytest.approx(1.0)

    def test_housekeeping_gene_has_low_tau(self) -> None:
        tau = tse.compute_tau(_synthetic_expression())
        # evenly-moderate [0,1,1,1] profile gives tau = 1/3
        assert tau.loc["og_housekeeping", "tau"] == pytest.approx(1.0 / 3.0)

    def test_silent_gene_is_nan(self) -> None:
        tau = tse.compute_tau(
            pd.DataFrame(
                {"a": [5.0, 0.0], "b": [0.0, 0.0]}, index=["g1", "silent"]
            )
        )
        assert tau.loc["g1", "tau"] == pytest.approx(1.0)
        assert np.isnan(tau.loc["silent", "tau"])

    def test_min_mean_expression_filter(self) -> None:
        expr = _synthetic_expression()
        tau = tse.compute_tau(expr, min_mean_expression=1.0)
        # og_housekeeping has mean 0.75 < 1.0 and is dropped.
        assert list(tau.index) == ["og_brain", "og_muscle", "og_gut"]
        assert tau.loc["og_brain", "tau"] == pytest.approx(1.0)

        silent = expr.copy()
        silent.loc["og_silent"] = 0.0
        tau_filtered = tse.compute_tau(silent, min_mean_expression=0.5)
        assert "og_silent" not in tau_filtered.index

    def test_requires_two_tissues(self) -> None:
        with pytest.raises(ValueError, match="two tissues"):
            tse.compute_tau(pd.DataFrame({"brain": [1.0]}, index=["g"]))

    def test_rejects_negative_expression(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            tse.compute_tau(pd.DataFrame({"a": [-1.0], "b": [1.0]}, index=["g"]))

    def test_scale_invariance(self) -> None:
        expr = _synthetic_expression()
        tau_a = tse.compute_tau(expr)
        tau_b = tse.compute_tau(expr * 1000.0)
        pd.testing.assert_frame_equal(tau_a, tau_b)


class TestAssignStates:
    def test_state_is_associated_tissue_when_tau_passes(self) -> None:
        expr = _synthetic_expression()
        tau = tse.compute_tau(expr)
        states = tse.assign_tissue_specificity_states(expr, tau, tau_cutoff=0.75)
        assert states.loc["og_brain", "tissue_specific_state"] == "brain"
        # tau = 1/3 does not pass the 0.75 cutoff.
        assert states.loc["og_housekeeping", "tissue_specific_state"] is None

    def test_below_cutoff_is_none(self) -> None:
        expr = _synthetic_expression()
        tau = tse.compute_tau(expr)
        states = tse.assign_tissue_specificity_states(expr, tau, tau_cutoff=0.95)
        assert states.loc["og_housekeeping", "tissue_specific_state"] is None
        assert states.loc["og_brain", "tissue_specific_state"] == "brain"

    def test_cutoff_boundary_inclusive(self) -> None:
        expr = pd.DataFrame(
            {"a": [2.0, 1.0], "b": [1.0, 1.0], "c": [1.0, 1.0]},
            index=["g", "other"],
        )
        tau = tse.compute_tau(expr)
        assert tau.loc["g", "tau"] == pytest.approx(0.5)
        states = tse.assign_tissue_specificity_states(expr, tau, tau_cutoff=0.5)
        assert states.loc["g", "tissue_specific_state"] == "a"


class TestParsimony:
    NEWICK = "((A:1,B:1):1,(C:1,D:1):1);"

    def test_uniform_state_zero_transitions(self) -> None:
        result = tse.count_parsimony_transitions(
            self.NEWICK, {"A": 1, "B": 1, "C": 1, "D": 1}
        )
        assert result["parsimony_score"] == 0
        assert result["root_state_candidates"] == [1]

    def test_single_tip_gain_is_one_transition(self) -> None:
        result = tse.count_parsimony_transitions(
            self.NEWICK, {"A": 1, "B": 0, "C": 0, "D": 0}
        )
        assert result["parsimony_score"] == 1
        assert result["root_state_candidates"] == [0]

    def test_two_independent_gains(self) -> None:
        result = tse.count_parsimony_transitions(
            self.NEWICK, {"A": 1, "B": 0, "C": 1, "D": 0}
        )
        # Either two gains on terminal branches or a root gain with two
        # losses; Fitch minimal score is 2 either way.
        assert result["parsimony_score"] == 2

    def test_clade_gain(self) -> None:
        result = tse.count_parsimony_transitions(
            self.NEWICK, {"A": 1, "B": 1, "C": 0, "D": 0}
        )
        assert result["parsimony_score"] == 1
        # Fitch root set is the union of alternatives: 0 or 1 both explain
        # the tips with a single change.
        assert set(result["root_state_candidates"]) == {0, 1}

    def test_missing_taxon_raises(self) -> None:
        with pytest.raises(ValueError, match="absent from tree"):
            tse.count_parsimony_transitions(
                self.NEWICK, {"A": 1, "B": 0, "C": 0, "E": 0}
            )

    def test_partial_tip_states_raise(self) -> None:
        with pytest.raises(ValueError, match="missing tip states"):
            tse.count_parsimony_transitions(self.NEWICK, {"A": 1, "B": 0})

    def test_five_taxon_asymmetric_tree(self) -> None:
        newick = "(((A:1,B:1):1,C:2):1,(D:1,E:1):2);"
        result = tse.count_parsimony_transitions(
            newick, {"A": 1, "B": 0, "C": 0, "D": 1, "E": 1}
        )
        assert result["parsimony_score"] == 2


class TestGainsLosses:
    NEWICK = "((apis:1,bombus:1):1,(ceratina:1,megachile:1):1);"

    def _states_frame(self, taus: dict[str, float]) -> pd.DataFrame:
        rows = [{"tau": tau, "species": species} for species, tau in taus.items()]
        frame = pd.DataFrame(rows, index=["og1"] * len(rows))
        frame.index.name = None
        return frame

    def test_full_observation_scored(self) -> None:
        states = self._states_frame(
            {"apis": 0.9, "bombus": 0.2, "ceratina": 0.2, "megachile": 0.2}
        )
        result = tse.count_tissue_specificity_gains_losses(self.NEWICK, states)
        assert result.loc["og1", "parsimony_score"] == 1
        assert result.loc["og1", "n_observed"] == 4
        assert result.loc["og1", "n_specific_species"] == 1
        assert bool(result.loc["og1", "descriptive"])

    def test_partial_observation_not_scored(self) -> None:
        states = self._states_frame({"apis": 0.9, "bombus": 0.2})
        result = tse.count_tissue_specificity_gains_losses(self.NEWICK, states)
        assert result.loc["og1", "parsimony_score"] is None
        assert result.loc["og1", "n_observed"] == 2

    def test_single_observation_not_scored(self) -> None:
        states = self._states_frame({"apis": 0.9})
        result = tse.count_tissue_specificity_gains_losses(self.NEWICK, states)
        assert result.loc["og1", "parsimony_score"] is None

    def test_cutoff_controls_state(self) -> None:
        states = self._states_frame(
            {"apis": 0.9, "bombus": 0.8, "ceratina": 0.2, "megachile": 0.2}
        )
        strict = tse.count_tissue_specificity_gains_losses(
            self.NEWICK, states, tau_cutoff=0.75
        )
        loose = tse.count_tissue_specificity_gains_losses(
            self.NEWICK, states, tau_cutoff=0.95
        )
        assert strict.loc["og1", "n_specific_species"] == 2
        assert loose.loc["og1", "n_specific_species"] == 0


class TestCoupling:
    def test_multicopy_enrichment_is_descriptive(self) -> None:
        rows = []
        taus = {"og_multi": {"A": 0.9, "B": 0.2, "C": 0.2, "D": 0.2},
                "og_single": {"A": 0.2, "B": 0.2, "C": 0.2, "D": 0.2}}
        index = []
        for og, per_species in taus.items():
            for species, tau in per_species.items():
                rows.append({"tau": tau, "species": species})
                index.append(og)
        gains = tse.count_tissue_specificity_gains_losses(
            "((A:1,B:1):1,(C:1,D:1):1);",
            pd.DataFrame(rows, index=index),
        )
        counts = {"og_multi": 3, "og_single": 1}
        table = tse.orthology_class_coupling_table(gains, counts)
        assert table.loc["multicopy", "n_with_transition"] == 1
        assert table.loc["single_copy", "n_with_transition"] == 0
        assert table.loc["multicopy", "fraction_with_transition"] == 1.0
        assert bool(table["descriptive"].all())
        assert "p_value" not in table.columns

    def test_missing_copy_count_defaults_single(self) -> None:
        gains = tse.count_tissue_specificity_gains_losses(
            "((A:1,B:1):1,(C:1,D:1):1);",
            pd.DataFrame(
                {
                    "tau": [0.9, 0.2, 0.2, 0.2],
                    "species": ["A", "B", "C", "D"],
                },
                index=["og_x"] * 4,
            ),
        )
        table = tse.orthology_class_coupling_table(gains, {})
        assert table.loc["single_copy", "n_orthogroups"] == 1
        assert table.loc["single_copy", "n_with_transition"] == 1


class TestTreeCoverage:
    NEWICK = "((apis:1,bombus:1):1,(ceratina:1,dummy:1):1);"

    def test_coverage_report(self) -> None:
        report = tse.validate_species_tree_coverage(
            self.NEWICK, ["apis", "bombus", "ceratina", "polistes"]
        )
        assert report["missing_species"] == ["polistes"]
        assert report["n_covered"] == 3
        assert report["n_requested"] == 4

    def test_full_coverage(self) -> None:
        report = tse.validate_species_tree_coverage(
            self.NEWICK, ["apis", "bombus", "ceratina", "dummy"]
        )
        assert report["missing_species"] == []
        assert report["n_covered"] == 4
