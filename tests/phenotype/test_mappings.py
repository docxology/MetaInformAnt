"""Tests for phenotype GWAS mapping constants and helpers."""

from __future__ import annotations

import pandas as pd
import pytest

from metainformant.phenotype.mappings import (
    BIOLOGICAL_GROUP_MAP,
    PHENOTYPE_LINK_MAP,
    STRAIN_ORDER,
    STRAIN_PALETTE,
    map_biological_groups,
    map_phenotype_links,
)


@pytest.fixture
def codes_df() -> pd.DataFrame:
    return pd.DataFrame({"biological_group_code": ["WORK", "ITW", "ITQ", "IV", "G"]})


class TestConstants:
    def test_group_and_link_maps_share_codes(self):
        assert set(BIOLOGICAL_GROUP_MAP) == set(PHENOTYPE_LINK_MAP)

    def test_strain_order_matches_palette(self):
        assert STRAIN_ORDER == list(STRAIN_PALETTE.keys())

    def test_phenotype_link_axis_values(self):
        assert set(PHENOTYPE_LINK_MAP.values()) == {"W", "Q"}
        # Workers: WORK and in vitro workers
        assert PHENOTYPE_LINK_MAP["WORK"] == "W"
        assert PHENOTYPE_LINK_MAP["ITW"] == "W"
        # Queens: in vitro queens, in vivo queens, grafted
        assert PHENOTYPE_LINK_MAP["ITQ"] == "Q"
        assert PHENOTYPE_LINK_MAP["IV"] == "Q"
        assert PHENOTYPE_LINK_MAP["G"] == "Q"


class TestMapBiologicalGroups:
    def test_adds_labels_column(self, codes_df):
        result = map_biological_groups(codes_df)
        assert result is codes_df  # documented in-place behaviour
        assert list(result["biological_group"]) == [
            "Worker",
            "In vitro worker",
            "In vitro queen",
            "In vivo (Queen)",
            "Grafted",
        ]

    def test_custom_code_column(self):
        df = pd.DataFrame({"bg": ["WORK"]})
        result = map_biological_groups(df, biological_group_code_column="bg")
        assert list(result["biological_group"]) == ["Worker"]

    def test_unknown_codes_map_to_nan(self):
        df = pd.DataFrame({"biological_group_code": ["WORK", "ZZZ"]})
        result = map_biological_groups(df)
        assert result["biological_group"].iloc[0] == "Worker"
        assert pd.isna(result["biological_group"].iloc[1])


class TestMapPhenotypeLinks:
    def test_adds_link_column(self, codes_df):
        result = map_phenotype_links(codes_df)
        assert result is codes_df
        assert list(result["phenotype_link"]) == ["W", "W", "Q", "Q", "Q"]

    def test_custom_code_column(self):
        df = pd.DataFrame({"code": ["ITQ"]})
        result = map_phenotype_links(df, biological_group_code_column="code")
        assert list(result["phenotype_link"]) == ["Q"]
