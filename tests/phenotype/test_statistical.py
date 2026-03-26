import pandas as pd
import pytest
from metainformant.phenotype.analysis.statistical import (
    calculate_summary_stats,
    perform_anova,
    perform_kruskal,
    perform_ttest,
    correlate_phenotypes
)

@pytest.fixture
def sample_df():
    data = {
        "group": ["A", "A", "A", "A", "B", "B", "B", "B"],
        "value": [10, 12, 11, 15, 20, 22, 21, 25],
        "value2": [100, 120, 110, 150, 200, 220, 210, 250]
    }
    return pd.DataFrame(data)

def test_calculate_summary_stats(sample_df):
    summary = calculate_summary_stats(sample_df, "value", "group")
    assert not summary.empty
    assert list(summary["group"]) == ["A", "B"]
    assert summary.iloc[0]["count"] == 4
    assert summary.iloc[0]["mean"] == 12.0
    assert summary.iloc[1]["mean"] == 22.0

def test_perform_anova(sample_df):
    res = perform_anova(sample_df, "value", "group")
    assert res.get("test") == "One-way ANOVA"
    assert "f_statistic" in res
    assert "p_value" in res
    assert res["groups"] == 2
    assert res["significant"] is True  # Values highly separated

def test_perform_kruskal(sample_df):
    res = perform_kruskal(sample_df, "value", "group")
    assert res.get("test") == "Kruskal-Wallis H-test"
    assert "h_statistic" in res
    assert "p_value" in res
    assert res["groups"] == 2
    assert res["significant"] is True

def test_perform_ttest(sample_df):
    res = perform_ttest(sample_df, "value", "group", "A", "B")
    assert res.get("test") == "Welch's t-test"
    assert res["group1"] == "A"
    assert res["group2"] == "B"
    assert "t_statistic" in res
    assert res["significant"] is True

def test_correlate_phenotypes(sample_df):
    corr = correlate_phenotypes(sample_df, ["value", "value2"])
    assert not corr.empty
    assert corr.shape == (2, 2)
    assert abs(corr.loc["value", "value2"] - 1.0) < 1e-6  # perfect correlation in dummy data

def test_invalid_inputs(sample_df):
    # Empty df
    res = perform_anova(pd.DataFrame(), "value", "group")
    assert "error" in res
    
    # Missing col
    res = perform_ttest(sample_df, "missing", "group", "A", "B")
    assert "error" in res
