import pandas as pd
import pytest

from metainformant.phenotype.analysis.statistical import (
    calculate_summary_stats,
    correlate_phenotypes,
    get_comprehensive_pairwise_ttests,
    perform_anova,
    perform_kruskal,
    perform_linear_regression,
    perform_multifactor_anova,
    perform_ttest,
)


@pytest.fixture
def sample_df():
    data = {
        "group": ["A", "A", "A", "A", "B", "B", "B", "B"],
        "value": [10, 12, 11, 15, 20, 22, 21, 25],
        "value2": [100, 120, 110, 150, 200, 220, 210, 250],
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


def test_perform_linear_regression(sample_df):
    # value2 is exactly 10 * value in the fixture data
    res = perform_linear_regression(sample_df, "value", "value2")
    assert "error" not in res
    assert res["test"] == "Linear Regression"
    assert res["slope"] == pytest.approx(10.0)
    assert res["r_squared"] == pytest.approx(1.0)
    assert res["n_obs"] == 8
    assert res["significant"] is True


def test_perform_linear_regression_errors(sample_df):
    res = perform_linear_regression(pd.DataFrame(), "value", "value2")
    assert "error" in res

    res = perform_linear_regression(sample_df, "missing", "value2")
    assert "error" in res

    # Fewer than 3 valid observations
    short = sample_df.iloc[:2]
    res = perform_linear_regression(short, "value", "value2")
    assert "error" in res


def test_perform_multifactor_anova():
    rng = pd.DataFrame(
        {
            "y": [1.0, 2.0, 3.0, 6.0, 7.0, 8.0, 2.0, 3.0, 4.0, 7.0, 8.0, 9.0],
            "f1": ["a", "a", "a", "a", "a", "a", "b", "b", "b", "b", "b", "b"],
            "f2": ["x", "x", "x", "y", "y", "y", "x", "x", "x", "y", "y", "y"],
        }
    )
    res = perform_multifactor_anova(rng, "y ~ C(f1) * C(f2)")
    assert "error" not in res
    assert set(res.keys()) == {"C(f1)", "C(f2)", "C(f1):C(f2)"}
    for row in res.values():
        assert "F" in row and "P_value" in row
        assert 0.0 <= row["P_value"] <= 1.0


def test_perform_multifactor_anova_unfittable_formula_returns_error():
    df = pd.DataFrame({"y": [1.0, 2.0, 3.0], "g": ["a", "b", "a"]})
    # Missing column makes the OLS fit itself fail; must return an error dict,
    # not raise NameError from the fallback path.
    res = perform_multifactor_anova(df, "y ~ C(not_a_column)")
    assert res == {"error": "Model could not be fit."}


def test_get_comprehensive_pairwise_ttests_sorted_and_complete(sample_df):
    results = get_comprehensive_pairwise_ttests(sample_df, "value", "group")
    assert len(results) == 1  # one pair: A vs B
    assert results[0]["group1"] == "A"
    assert results[0]["group2"] == "B"

    df3 = pd.concat([sample_df, sample_df.assign(group="C", value=sample_df["value"] + 40)], ignore_index=True)
    results = get_comprehensive_pairwise_ttests(df3, "value", "group")
    assert len(results) == 3
    p_values = [r["p_value"] for r in results]
    assert p_values == sorted(p_values)


def test_get_comprehensive_pairwise_ttests_skips_empty_groups(sample_df):
    # A group whose values are all NaN yields no valid t-test and is skipped
    df = pd.concat([sample_df, sample_df.assign(group="C", value=float("nan"))], ignore_index=True)
    results = get_comprehensive_pairwise_ttests(df, "value", "group")
    assert len(results) == 1
    assert {results[0]["group1"], results[0]["group2"]} == {"A", "B"}
