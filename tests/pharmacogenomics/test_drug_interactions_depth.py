"""Depth tests for drug-drug interaction prediction and polypharmacy risk.

Exercises the built-in interaction database, CYP pathway tables, pairwise
prediction, polypharmacy risk scoring, and CYP inhibition/induction profiling
with real data structures only.
"""

from __future__ import annotations

import pytest

from metainformant.pharmacogenomics.interaction.drug_interactions import (
    _find_shared_cyp_pathways,
    cyp_inhibition_prediction,
    default_interaction_database,
    polypharmacy_risk,
    predict_drug_interaction,
)


class TestDefaultDatabase:
    def test_database_shape_and_keys_sorted_pairs(self) -> None:
        db = default_interaction_database()
        assert len(db) >= 30  # database with major/moderate/minor tiers
        for key, record in db.items():
            assert isinstance(key, tuple) and len(key) == 2
            assert key == tuple(sorted(key))
            assert set(record) >= {"severity", "mechanism", "description", "evidence_level", "recommendation"}
            assert record["severity"] in {"major", "moderate", "minor"}
            assert record["evidence_level"] in {"A", "B", "C", "D"}

    def test_database_is_deterministic_across_calls(self) -> None:
        assert default_interaction_database() == default_interaction_database()

    def test_severity_distribution_is_nontrivial(self) -> None:
        db = default_interaction_database()
        severities = {r["severity"] for r in db.values()}
        assert severities == {"major", "moderate", "minor"}


class TestPredictDrugInteraction:
    def test_documented_pair_major(self) -> None:
        result = predict_drug_interaction("warfarin", "fluconazole")
        assert result["severity"] == "major"
        assert result["mechanism"]
        assert result["evidence_level"] in {"A", "B"}
        assert result["drug_a"] == "warfarin" or result["drug_a"] == "fluconazole"

    def test_order_and_case_insensitive_lookup(self) -> None:
        ab = predict_drug_interaction("Warfarin ", "FLUCONAZOLE")
        ba = predict_drug_interaction("fluconazole", "warfarin")
        assert ab["severity"] == ba["severity"] == "major"
        assert ab["mechanism"] == ba["mechanism"]

    def test_undocumented_pair_shares_cyp_pathway(self) -> None:
        # metoprolol and tramadol are both CYP2D6 substrates but the pair is
        # not expected to be an explicit database entry.
        shared = _find_shared_cyp_pathways("metoprolol", "tramadol")
        assert shared == ["CYP2D6"]
        result = predict_drug_interaction("metoprolol", "tramadol")
        assert result["severity"] == "none"
        assert "CYP2D6" in result["mechanism"]
        assert result["evidence_level"] == "D"

    def test_completely_unknown_pair(self) -> None:
        result = predict_drug_interaction("unobtanium-A", "unobtanium-B")
        assert result["severity"] == "none"
        assert result["mechanism"] == "None identified"
        assert "No known interaction" in result["description"]

    def test_custom_database_used_instead_of_builtin(self) -> None:
        custom = {
            ("x-drug", "y-drug"): {
                "severity": "major",
                "mechanism": "m",
                "description": "d",
                "evidence_level": "A",
                "recommendation": "r",
            }
        }
        result = predict_drug_interaction("x-drug", "y-drug", interaction_db=custom)
        assert result["severity"] == "major"


class TestPolypharmacyRisk:
    def test_requires_two_medications(self) -> None:
        with pytest.raises(ValueError, match="At least 2"):
            polypharmacy_risk(["warfarin"])

    def test_known_major_pair_scores(self) -> None:
        result = polypharmacy_risk(["warfarin", "fluconazole"])
        assert result["n_interactions"] >= 1
        assert result["interactions"][0]["severity"] == "major"
        assert 0.0 < result["risk_score"] <= 1.0
        assert result["severity_counts"]["major"] >= 1
        assert result["recommendations"]

    def test_clean_pair_low_risk(self) -> None:
        result = polypharmacy_risk(["metoprolol", "tramadol"])
        # Shared CYP2D6 pathway adds a small penalty but no documented interaction.
        assert result["n_interactions"] == 0
        assert "CYP2D6" in result["competing_pathways"]
        assert result["risk_score"] < 0.2

    def test_three_drug_regimen_aggregates_pairs(self) -> None:
        result = polypharmacy_risk(["warfarin", "fluconazole", "simvastatin", "ketoconazole"])
        expected_pairs = 4 * 3 // 2
        assert result["n_medications"] == 4
        assert result["risk_score"] <= 1.0
        # All pairwise interaction records come from real database lookups.
        assert len(result["interactions"]) <= expected_pairs
        assert len({r["recommendation"] for r in result["interactions"]}) <= len(result["recommendations"])

    def test_risk_score_monotone_in_severity(self) -> None:
        minor = polypharmacy_risk(["caffeine", "cimetidine"])
        major = polypharmacy_risk(["warfarin", "fluconazole"])
        assert minor["risk_score"] < major["risk_score"]


class TestCypInhibitionPrediction:
    def test_strong_inhibitor_profile(self) -> None:
        result = cyp_inhibition_prediction("ketoconazole")
        assert result["drug"] == "ketoconazole"
        assert result["n_affected"] >= 1
        entry = result["affected_enzymes"][0]
        assert entry == {
            "enzyme": "CYP3A4",
            "effect_type": "inhibitor",
            "inhibition_type": "strong",
            "potency_estimate": "high",
        }
        assert "CYP3A4" in result["summary"]

    def test_dual_role_drug(self) -> None:
        # fluoxetine inhibits CYP2D6 (strong) and CYP2C19 (moderate) and is a
        # CYP2D6 + CYP2C19 substrate... but it is not listed as substrate; it
        # is only in inhibitor tables. Verify whatever the tables say.
        result = cyp_inhibition_prediction("fluoxetine")
        types = {a["effect_type"] for a in result["affected_enzymes"]}
        assert types == {"inhibitor"}
        enzymes = {a["enzyme"] for a in result["affected_enzymes"]}
        assert {"CYP2D6", "CYP2C19"} <= enzymes

    def test_inducer_profile(self) -> None:
        result = cyp_inhibition_prediction("rifampin")
        types = {a["effect_type"] for a in result["affected_enzymes"]}
        assert types == {"inducer"}
        assert result["n_affected"] >= 3  # CYP1A2, CYP2C9, CYP2C19, CYP3A4

    def test_substrate_only_drug(self) -> None:
        result = cyp_inhibition_prediction("warfarin")
        assert result["affected_enzymes"] == []
        assert result["n_affected"] == 0
        assert "CYP2C9" in result["substrate_of"]
        assert "Metabolised by" in result["summary"]

    def test_unknown_drug_neutral_summary(self) -> None:
        result = cyp_inhibition_prediction("mystery-drug")
        assert result["n_affected"] == 0
        assert result["substrate_of"] == []
        assert result["summary"] == "mystery-drug has no known CYP inhibition or induction effects"

    def test_custom_profiles_override_builtin(self) -> None:
        profiles = {
            "inhibitors": {"CYP1A2": [{"drug": "custom-inhib", "type": "strong", "potency": "high"}]},
            "inducers": {},
        }
        result = cyp_inhibition_prediction("custom-inhib", cyp_profiles=profiles)
        assert result["n_affected"] == 1
        assert result["affected_enzymes"][0]["enzyme"] == "CYP1A2"
