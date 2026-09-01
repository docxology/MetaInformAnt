"""Tests for find_differential_pathways (metagenomics functional pathways).

Completes coverage of the pathways module public API. Uses real synthetic
annotation data and real implementations only (real-implementation policy).
"""

from __future__ import annotations

import pytest

from metainformant.metagenomics.functional.pathways import (
    _BUILTIN_PATHWAYS,
    PathwayDefinition,
    find_differential_pathways,
)

_GLYCOLYSIS_KOS = [
    "K00844",
    "K00845",
    "K01810",
    "K00850",
    "K01623",
    "K01624",
    "K01803",
    "K00134",
    "K00927",
    "K01689",
    "K00873",
]

_TCA_KOS = [
    "K01647",
    "K01681",
    "K00031",
    "K00030",
    "K00164",
    "K00658",
    "K01902",
    "K01903",
    "K00234",
    "K00235",
    "K01676",
    "K00024",
]


def _annotations_with(kos: list[str], fillers: int = 0) -> dict[str, list[str]]:
    """One gene per KO plus optional filler genes with unrelated annotations."""
    ann = {f"gene_{ko}": [ko] for ko in kos}
    for i in range(fillers):
        ann[f"filler_{i}"] = [f"K9999{i}"]
    return ann


class TestFindDifferentialPathways:
    def test_group1_enriched_glycolysis(self) -> None:
        group1 = {f"s{i}": _annotations_with(_GLYCOLYSIS_KOS) for i in range(3)}
        group2 = {f"s{i}": _annotations_with(_TCA_KOS) for i in range(3)}

        results = find_differential_pathways(group1, group2, min_diff=0.2)

        by_id = {r["pathway_id"]: r for r in results}
        assert "map00010" in by_id
        row = by_id["map00010"]
        assert row["group1_mean"] == 1.0
        assert row["group2_mean"] == 0.0
        assert row["difference"] == 1.0
        assert row["direction"] == "group1_enriched"
        assert row["pathway_name"] == "Glycolysis / Gluconeogenesis"

    def test_min_diff_threshold_filters(self) -> None:
        group1 = {f"s{i}": _annotations_with(_GLYCOLYSIS_KOS) for i in range(2)}
        group2 = {f"s{i}": _annotations_with(_TCA_KOS) for i in range(2)}
        # A threshold above the max possible |diff| (1.0) yields nothing.
        assert find_differential_pathways(group1, group2, min_diff=1.5) == []

    def test_custom_definitions_use_given_name(self) -> None:
        custom = {
            "pw_x": PathwayDefinition(
                pathway_id="pw_x",
                name="Custom Pathway X",
                required_kos=["K00001", "K00002", "K00003", "K00004"],
            )
        }
        group1 = {
            "a": {"g1": ["K00001", "K00002", "K00003", "K00004"]},
            "b": {"g2": ["K00001", "K00002", "K00003", "K00004"]},
        }
        group2 = {"c": {"g3": ["K00001"]}, "d": {"g4": ["K00001"]}}

        results = find_differential_pathways(group1, group2, pathway_definitions=custom, min_diff=0.2)

        assert len(results) == 1
        row = results[0]
        assert row["pathway_id"] == "pw_x"
        assert row["pathway_name"] == "Custom Pathway X"
        assert row["group1_mean"] == pytest.approx(1.0)
        assert row["group2_mean"] == pytest.approx(0.25)  # 1 of 4 required KOs

    def test_intermediate_completeness_means(self) -> None:
        # Group 1 has full glycolysis; group 2 has exactly half the required KOs.
        half = _GLYCOLYSIS_KOS[:5]  # 5 of 11 required KOs
        group1 = {"s1": _annotations_with(_GLYCOLYSIS_KOS)}
        group2 = {"s2": _annotations_with(half)}

        results = find_differential_pathways(group1, group2, min_diff=0.2)
        by_id = {r["pathway_id"]: r for r in results}
        row = by_id["map00010"]
        assert row["group2_mean"] == pytest.approx(5 / 11, abs=1e-4)  # module rounds to 4 dp
        assert row["difference"] == pytest.approx(1.0 - 5 / 11, abs=1e-4)

    def test_builtin_names_resolved_without_custom_defs(self) -> None:
        group1 = {"s1": _annotations_with(_TCA_KOS)}
        group2 = {"s2": {"x": ["K99999"]}}
        results = find_differential_pathways(group1, group2, min_diff=0.2)
        by_id = {r["pathway_id"]: r for r in results}
        assert by_id["map00020"]["pathway_name"] == "Citrate cycle (TCA cycle)"
        assert "map00020" in _BUILTIN_PATHWAYS
