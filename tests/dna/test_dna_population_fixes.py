"""Regression tests for dna.population fixes (LD pairing, Fay & Wu H, Fst filtering, plotting)."""

# ruff: noqa: E402 - Agg backend must be forced before pyplot-backed modules are imported.

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from metainformant.dna.population import (  # noqa: E402
    analysis,
    core,
    visualization_core as vizcore,
    visualization_stats as vizstats,
)


class TestLinkageDisequilibrium:
    """linkage_disequilibrium must pair haplotypes from the same sequence."""

    def test_hand_computed_d_and_dprime(self):
        # Haplotypes: (A,T) x3 and (G,C) x1 -> complete linkage between the
        # minor alleles.  Using the implementation's convention (majority base
        # is ancestral): pA = 3/4, pB = 3/4, pAB = 3/4.
        seqs = ["AT", "AT", "AT", "GC"]
        d = core.linkage_disequilibrium(seqs, 0, 1)
        assert d == pytest.approx(0.75 - 0.75 * 0.75)  # 0.1875

        # D' = D / D_max with D_max = min(pA*(1-pB), (1-pA)*pB) = 0.1875 -> 1.0
        p_a, p_b = 0.75, 0.75
        d_max = min(p_a * (1 - p_b), (1 - p_a) * p_b)
        assert d / d_max == pytest.approx(1.0)

    def test_sequence_ambiguous_at_one_position_excluded_entirely(self):
        # "AN" is valid at position 0 but ambiguous at position 1: it must not
        # contribute an allele to either column, and must not shift pairing.
        clean = ["AT", "AT", "AT", "GC"]
        noisy = clean + ["AN"]
        expected = 0.75 - 0.75 * 0.75
        assert core.linkage_disequilibrium(noisy, 0, 1) == pytest.approx(expected)

    def test_pairing_shift_from_mixed_ambiguity_is_gone(self):
        # "NT" is valid only at position 1 and "AN" only at position 0.  With
        # the old independent filtering both columns ended up with n=3 and the
        # haplotype pairing shifted; now only fully valid sequences pair.
        noisy = ["AT", "NT", "AN", "GC"]
        kept_only = ["AT", "GC"]
        assert core.linkage_disequilibrium(noisy, 0, 1) == pytest.approx(core.linkage_disequilibrium(kept_only, 0, 1))


class TestFayWuH:

    def test_hand_computed_per_site_h(self):
        seqs = ["AAAA", "AAAC", "AGGC", "TGGC"]
        outgroup = "AAAA"
        # Derived-allele counts per site (outgroup = AAAA is ancestral):
        # site 0: i=1 -> 2*1/12; sites 1, 2: i=2 -> 2*4/12 each; site 3: i=3 -> 2*9/12.
        theta_h_raw = 2 * 1 / 12 + 2 * 4 / 12 + 2 * 4 / 12 + 2 * 9 / 12  # 3.0
        theta_h_per_site = theta_h_raw / 4
        # Pairwise differences: 1,3,4,2,3,1 of 4 sites -> pi = (14/4)/6
        pi = ((1 + 3 + 4 + 2 + 3 + 1) / 4) / 6
        assert core.fay_wu_h_from_sequences(seqs, outgroup) == pytest.approx(pi - theta_h_per_site)

    def test_monomorphic_flanks_dilute_h_by_exact_length_factor(self):
        # Per-site quantities scale exactly by L/(L+f) when f monomorphic
        # sites are appended, so H must scale by the same factor.
        seqs = ["AAAA", "AAAC", "AGGC", "TGGC"]
        outgroup = "AAAA"
        flank = "GG"
        h_short = core.fay_wu_h_from_sequences(seqs, outgroup)
        h_long = core.fay_wu_h_from_sequences([s + flank for s in seqs], outgroup + flank)
        assert h_long == pytest.approx(h_short * 4 / 6)


class TestCalculateFst:
    """calculate_fst must ignore gaps/ambiguous bases, not count them as alleles."""

    def test_gap_monomorphic_site_does_not_change_fst(self):
        # Site 0 is polymorphic (A/A vs G/A); site 1 is A,-,A,A: a gap next to
        # a monomorphic base must not make the site polymorphic.
        with_gap = (["AA", "A-"], ["GA", "AA"])
        without_site = (["A", "A"], ["G", "A"])
        assert analysis.calculate_fst(*with_gap) == pytest.approx(analysis.calculate_fst(*without_site))

    def test_hand_computed_fst_value(self):
        # Site 0: pop1 = A,A and pop2 = G,A. ht = 0.375, hs = 0.25 -> 1/3.
        assert analysis.calculate_fst(["A", "A"], ["G", "A"]) == pytest.approx(1.0 / 3.0)


class TestEstimatePopulationSize:
    """Ne must follow ne = theta_w / (4 * mu * L) with theta_w = S/a_n."""

    def test_hand_computed_ne(self):
        seqs = ["ATCG", "ATCG", "GCTA"]  # S = 4, n = 3, a_n = 1.5
        sizes = analysis.estimate_population_size(seqs)
        theta_w = 4 / 1.5
        assert sizes["theta_watterson"] == pytest.approx(theta_w)
        assert sizes["ne_estimate"] == pytest.approx(theta_w / (4 * 1e-8 * 4))


class TestInterpretNeutralityResults:
    """Demographic labels must match the sign of Tajima's D."""

    def test_positive_d_labels(self):
        result = analysis.interpret_neutrality_results({"tajima_d": 1.0})
        assert result["tajima_d"] == "balancing_selection_or_population_bottleneck"

    def test_negative_d_labels(self):
        result = analysis.interpret_neutrality_results({"tajima_d": -1.0})
        assert result["tajima_d"] == "positive_selection_or_population_expansion"

    def test_missing_value_labels(self):
        result = analysis.interpret_neutrality_results({"tajima_d": None})
        assert result["tajima_d"] == "calculation_failed"


def _scenario_fst_matrix():
    return np.array([[0.0, 0.1], [0.1, 0.0]])


def _scenario_pca_coords():
    return np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.5], [0.5, 2.0]])


PLOT_SCENARIOS = {
    "plot_fst_matrix": lambda: vizcore.plot_fst_matrix(_scenario_fst_matrix()),
    "plot_tajima_d_distribution": lambda: vizcore.plot_tajima_d_distribution([0.1, -0.2, 0.3, -0.1]),
    "plot_selection_statistics": lambda: vizcore.plot_selection_statistics({"tajima_d": [0.1, -0.2, 0.3]}),
    "plot_population_diversity": lambda: vizcore.plot_population_diversity({"pop1": 0.01, "pop2": 0.02}),
    "plot_ld_decay": lambda: vizcore.plot_ld_decay([(1, 0.8), (2, 0.6), (3, 0.4), (4, 0.3)]),
    "plot_population_structure": lambda: vizcore.plot_population_structure(_scenario_pca_coords()),
    "plot_demographic_history": lambda: vizcore.plot_demographic_history([1000.0, 1200.0], [0, 1000]),
    "create_population_summary_plot": lambda: vizcore.create_population_summary_plot(
        {"pop1": {"diversity": 0.01}, "pop2": {"diversity": 0.02}}
    ),
    "plot_mutation_spectrum": lambda: vizcore.plot_mutation_spectrum({"A>T": 10, "C>G": 5}),
    "plot_allele_frequency_spectrum": lambda: vizcore.plot_allele_frequency_spectrum([0.1, 0.2, 0.3, 0.8]),
    "plot_bootstrap_distribution": lambda: vizcore.plot_bootstrap_distribution([1.0, 1.1, 0.9, 1.2]),
    "plot_demographic_comparison": lambda: vizstats.plot_demographic_comparison(
        {"era1": {"estimated_ne": 1000, "observed_diversity": 0.5}}
    ),
    "plot_fst_comparison": lambda: vizstats.plot_fst_comparison({"locus1": 0.02, "locus2": 0.12}),
    "plot_hardy_weinberg_test": lambda: vizstats.plot_hardy_weinberg_test([{"locus": "L1", "p_value": 0.5}]),
    "plot_heterozygosity_distribution": lambda: vizstats.plot_heterozygosity_distribution([0.1, 0.2, 0.3]),
    "plot_kinship_matrix": lambda: vizstats.plot_kinship_matrix(np.array([[1.0, 0.5], [0.5, 1.0]])),
    "plot_linkage_disequilibrium_decay": lambda: vizstats.plot_linkage_disequilibrium_decay([(1, 0.8), (2, 0.6)]),
    "plot_neutrality_test_suite": lambda: vizstats.plot_neutrality_test_suite({"tajima_d": -1.0, "fu_li_d": 0.5}),
    "plot_neutrality_test_summary": lambda: vizstats.plot_neutrality_test_summary(
        {"tajima_d": {"statistic": -1.0, "p_value": 0.2}, "fu_li_d": {"statistic": 0.5}}
    ),
    "plot_outlier_detection": lambda: vizstats.plot_outlier_detection([0.1, 0.2, 9.9], outliers=[2]),
    "plot_permutation_test": lambda: vizstats.plot_permutation_test([0.1, 0.2, 0.3], observed_value=0.15, p_value=0.04),
    "plot_pi_vs_theta": lambda: vizstats.plot_pi_vs_theta([0.1, 0.2], [0.15, 0.25]),
    "plot_site_frequency_spectrum": lambda: vizstats.plot_site_frequency_spectrum([5, 3, 1]),
    "plot_statistic_correlation_matrix": lambda: vizstats.plot_statistic_correlation_matrix(
        {"a": [1.0, 2.0, 3.0], "b": [2.0, 4.0, 6.5]}
    ),
    "plot_statistic_distribution": lambda: vizstats.plot_statistic_distribution([1.0, 2.0, 3.0]),
    "plot_summary_statistics_grid": lambda: vizstats.plot_summary_statistics_grid(
        {"pop1": {"diversity": 0.01, "fst": 0.05}, "pop2": {"diversity": 0.02, "fst": 0.03}}
    ),
    "plot_tajimas_d_comparison": lambda: vizstats.plot_tajimas_d_comparison({"pop1": -1.0, "pop2": 2.5}),
}

FORMAT_STRING_PLOTTERS = (
    "plot_population_diversity",
    "plot_ld_decay",
    "plot_demographic_history",
    "create_population_summary_plot",
)

LITERAL_FORMAT_STRINGS = {".4f", ".3f", ".0f"}


@pytest.mark.parametrize("name", sorted(PLOT_SCENARIOS))
def test_plotter_returns_closed_figure(name: str):
    """Every population plotter returns a Figure and leaks nothing into pyplot."""
    fig = PLOT_SCENARIOS[name]()
    assert fig is not None
    assert len(plt.get_fignums()) == 0
    plt.close(fig)


@pytest.mark.parametrize("name", FORMAT_STRING_PLOTTERS)
def test_format_string_plotters_render_values_not_literals(name: str):
    """The fixed plotters must render formatted values, never the literal mask."""
    fig = PLOT_SCENARIOS[name]()
    assert fig is not None
    rendered = set()
    for ax in fig.axes:
        rendered.update(text.get_text() for text in ax.texts)
        if ax.get_legend() is not None:
            rendered.update(text.get_text() for text in ax.get_legend().get_texts())
    assert not rendered & LITERAL_FORMAT_STRINGS
    if name == "plot_ld_decay":
        # The formatted decay-rate label lives in the legend.
        assert any(text.startswith("Exponential fit (b=") for text in rendered)
    else:
        assert any(text for text in rendered)  # values were actually drawn
    plt.close(fig)


def test_plot_pca_results_draws_real_scatter():
    pca_result = {
        "status": "success",
        "pcs": [[0.0, 0.0], [1.0, 1.0], [2.0, 0.5], [0.5, 2.0]],
        "explained_variance_ratio": [0.6, 0.3],
        "labels": [0, 0, 1, 1],
    }
    fig = vizstats.plot_pca_results(pca_result)
    assert fig is not None
    assert len(fig.axes) == 1
    scatters = [collection for collection in fig.axes[0].collections if hasattr(collection, "get_offsets")]
    assert len(scatters) >= 1
    point_counts = [len(collection.get_offsets()) for collection in scatters]
    assert 4 in point_counts
    assert "PC1 (60.0% variance)" == fig.axes[0].get_xlabel()
    assert len(plt.get_fignums()) == 0
    plt.close(fig)


def test_plot_pca_results_returns_none_without_coordinates():
    assert vizstats.plot_pca_results({"status": "success"}) is None


def test_plot_pi_vs_theta_returns_none_on_empty_input():
    assert vizstats.plot_pi_vs_theta([], []) is None
    assert vizstats.plot_pi_vs_theta([0.1], []) is None
    assert len(plt.get_fignums()) == 0
