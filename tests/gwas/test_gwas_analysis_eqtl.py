"""Regression tests: ``run_eqtl_analysis`` applies its ``cis_window`` to variant and gene positions.

Defect: ``cis_window`` was documented as "not used in simple linear test, placeholder": the function tested every
variant against every transcript whatever window the caller passed, and had no way to be told where anything is.
It now accepts ``variant_positions`` and ``gene_positions`` (the tables ``cis_eqtl_scan`` uses) and restricts the
scan to variants on the transcript's chromosome within ``cis_window`` bp of its TSS. A ``cis_window`` without
positions is an error instead of being ignored.

Oracles: the expected cis pairs are listed by hand from the coordinates in ``VARIANTS`` and ``GENES`` (including a
pair that is close in base pairs but on another chromosome, and pairs exactly on the window edge); every reported
effect is compared with ``scipy.stats.linregress`` on the same two columns.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from metainformant.gwas.analysis.eqtl import run_eqtl_analysis

N_SAMPLES = 12
ONE_MEGABASE = 1_000_000  # the documented default window when positions are given and cis_window is omitted

# variant id -> (chromosome, position)
VARIANTS = {
    "V1": ("chr1", 1_000),
    "V2": ("chr1", 5_000),
    "V3": ("chr1", 20_000),
    "V4": ("chr1", 90_000),
    "V5": ("chr2", 1_500),
    "V6": ("chr2", 200_000),
}
# transcript id -> (chromosome, TSS)
GENES = {
    "T1": ("chr1", 4_000),
    "T2": ("chr2", 2_000),
    "T3": ("chr1", 95_000),
}


def _positions() -> tuple[pd.DataFrame, pd.DataFrame]:
    variant_positions = pd.DataFrame(
        [{"variant_id": vid, "chrom": chrom, "position": pos} for vid, (chrom, pos) in VARIANTS.items()]
    )
    gene_positions = pd.DataFrame(
        [{"gene_id": gid, "chrom": chrom, "tss_position": tss} for gid, (chrom, tss) in GENES.items()]
    )
    return variant_positions, gene_positions


@pytest.fixture(scope="module")
def data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Genotypes (samples x variants) and expression (samples x transcripts) with some real signal."""
    rng = np.random.default_rng(20261027)
    samples = [f"S{i}" for i in range(N_SAMPLES)]
    genotypes = pd.DataFrame({vid: rng.integers(0, 3, N_SAMPLES) for vid in VARIANTS}, index=samples)
    assert (genotypes.nunique() > 1).all()  # no constant column: every pair is estimable
    expression = pd.DataFrame(
        {
            "T1": 0.9 * genotypes["V2"] + rng.normal(0, 0.5, N_SAMPLES),
            "T2": -0.7 * genotypes["V5"] + rng.normal(0, 0.5, N_SAMPLES),
            "T3": 0.4 * genotypes["V4"] + rng.normal(0, 0.5, N_SAMPLES),
        },
        index=samples,
    )
    return genotypes, expression


def _pairs(result: pd.DataFrame) -> set[tuple[str, str]]:
    return set(zip(result["transcript"], result["variant"]))


def _expected_pairs(window: int) -> set[tuple[str, str]]:
    """Cis pairs by the plainest possible rule: same chromosome and |position - TSS| <= window."""
    return {
        (gene, variant)
        for gene, (gene_chrom, tss) in GENES.items()
        for variant, (variant_chrom, position) in VARIANTS.items()
        if gene_chrom == variant_chrom and abs(position - tss) <= window
    }


def test_the_hand_listed_cis_pairs_match_the_plain_rule() -> None:
    """Guard the oracle: the by-hand lists below are what the plain rule gives."""
    # T1 (chr1:4000) reaches V1 (3000 away, on the edge) and V2 (1000); V5 is 2500 away in bp but on chr2.
    # T2 (chr2:2000) reaches V5 (500); V1 is 1000 away in bp but on chr1. T3 (chr1:95000) reaches nothing.
    assert _expected_pairs(3_000) == {("T1", "V1"), ("T1", "V2"), ("T2", "V5")}
    # Widening to 5000 brings in T3-V4, exactly on the edge (|90000 - 95000| = 5000).
    assert _expected_pairs(5_000) == {("T1", "V1"), ("T1", "V2"), ("T2", "V5"), ("T3", "V4")}
    assert _expected_pairs(0) == set()


@pytest.mark.parametrize("window", [3_000, 5_000, 100_000, ONE_MEGABASE])
def test_only_pairs_inside_the_window_are_tested(data: tuple[pd.DataFrame, pd.DataFrame], window: int) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    result = run_eqtl_analysis(
        genotypes,
        expression,
        cis_window=window,
        variant_positions=variant_positions,
        gene_positions=gene_positions,
    )

    assert _pairs(result) == _expected_pairs(window)
    assert len(result) == len(_pairs(result))  # each pair is tested once


def test_the_default_window_is_one_megabase_when_positions_are_given(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    result = run_eqtl_analysis(
        genotypes, expression, variant_positions=variant_positions, gene_positions=gene_positions
    )

    assert _pairs(result) == _expected_pairs(ONE_MEGABASE)
    assert len(_expected_pairs(ONE_MEGABASE)) == 10  # 4 + 2 + 4: every same-chromosome pair, none across chromosomes


def test_the_window_edge_is_inclusive_and_distances_are_signed(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    result = run_eqtl_analysis(
        genotypes, expression, cis_window=3_000, variant_positions=variant_positions, gene_positions=gene_positions
    )

    distance = {(row.transcript, row.variant): row.distance for row in result.itertuples()}
    assert distance == {("T1", "V1"): -3_000, ("T1", "V2"): 1_000, ("T2", "V5"): -500}


def test_cross_chromosome_pairs_are_never_tested_however_close_in_base_pairs(
    data: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    result = run_eqtl_analysis(
        genotypes, expression, cis_window=10**9, variant_positions=variant_positions, gene_positions=gene_positions
    )

    assert ("T1", "V5") not in _pairs(result)  # V5 is 2500 bp from T1's TSS numerically, but on chr2
    assert ("T2", "V1") not in _pairs(result)
    assert all(VARIANTS[v][0] == GENES[t][0] for t, v in _pairs(result))


def test_chromosome_labels_are_compared_as_text(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    """Integer chromosome codes in one table and string codes in the other still match."""
    genotypes, expression = data
    variant_positions, gene_positions = _positions()
    variant_positions["chrom"] = variant_positions["chrom"].str.replace("chr", "").astype(int)
    gene_positions["chrom"] = gene_positions["chrom"].str.replace("chr", "")

    result = run_eqtl_analysis(
        genotypes, expression, cis_window=3_000, variant_positions=variant_positions, gene_positions=gene_positions
    )

    assert _pairs(result) == _expected_pairs(3_000)


def test_reported_effects_match_an_independent_regression(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    result = run_eqtl_analysis(
        genotypes, expression, cis_window=5_000, variant_positions=variant_positions, gene_positions=gene_positions
    )

    assert len(result) == 4
    for row in result.itertuples():
        oracle = stats.linregress(genotypes[row.variant], expression[row.transcript])
        assert row.beta == pytest.approx(oracle.slope, rel=1e-9)
        assert row.p_value == pytest.approx(oracle.pvalue, rel=1e-6)
    # The planted effects are the strongest signals of their transcripts.
    by_pair = {(row.transcript, row.variant): row.p_value for row in result.itertuples()}
    assert by_pair[("T1", "V2")] < by_pair[("T1", "V1")]


def test_a_cis_scan_reports_the_same_statistics_as_the_all_vs_all_scan(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    everything = run_eqtl_analysis(genotypes, expression)
    cis = run_eqtl_analysis(
        genotypes, expression, cis_window=5_000, variant_positions=variant_positions, gene_positions=gene_positions
    )

    assert len(everything) == len(VARIANTS) * len(GENES)
    assert "distance" not in everything.columns  # nothing to measure a distance against without positions
    by_pair = everything.set_index(["transcript", "variant"])
    for row in cis.itertuples():
        reference = by_pair.loc[(row.transcript, row.variant)]
        assert row.beta == pytest.approx(reference["beta"], rel=1e-12)
        assert row.p_value == pytest.approx(reference["p_value"], rel=1e-12)


def test_without_positions_every_variant_is_still_tested_against_every_transcript(
    data: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    genotypes, expression = data

    result = run_eqtl_analysis(genotypes, expression)

    assert _pairs(result) == {(g, v) for g in GENES for v in VARIANTS}


def test_max_transcripts_limits_a_cis_scan_too(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    result = run_eqtl_analysis(
        genotypes,
        expression,
        cis_window=5_000,
        max_transcripts=2,
        variant_positions=variant_positions,
        gene_positions=gene_positions,
    )

    assert _pairs(result) == {("T1", "V1"), ("T1", "V2"), ("T2", "V5")}  # T3 is the third transcript


def test_variants_and_transcripts_without_a_position_are_not_tested(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()
    variant_positions = variant_positions[variant_positions["variant_id"] != "V2"]  # V2 cannot be placed
    gene_positions = gene_positions[gene_positions["gene_id"] != "T2"]  # T2 cannot be placed

    result = run_eqtl_analysis(
        genotypes, expression, cis_window=5_000, variant_positions=variant_positions, gene_positions=gene_positions
    )

    assert _pairs(result) == {("T1", "V1"), ("T3", "V4")}


def test_a_window_without_positions_is_an_error_not_a_silent_no_op(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    """The defect itself: cis_window used to be accepted and ignored, returning the full all-vs-all table."""
    genotypes, expression = data

    with pytest.raises(ValueError, match="cis_window needs variant_positions and gene_positions"):
        run_eqtl_analysis(genotypes, expression, cis_window=3_000)


def test_the_two_position_tables_must_be_given_together(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    with pytest.raises(ValueError, match="must be given together"):
        run_eqtl_analysis(genotypes, expression, variant_positions=variant_positions)
    with pytest.raises(ValueError, match="must be given together"):
        run_eqtl_analysis(genotypes, expression, cis_window=3_000, gene_positions=gene_positions)


def test_a_negative_window_is_rejected(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    with pytest.raises(ValueError, match="non-negative"):
        run_eqtl_analysis(
            genotypes, expression, cis_window=-1, variant_positions=variant_positions, gene_positions=gene_positions
        )


def test_position_tables_need_their_columns_and_unique_identifiers(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    with pytest.raises(ValueError, match=r"gene_positions is missing required column\(s\) \['tss_position'\]"):
        run_eqtl_analysis(
            genotypes,
            expression,
            variant_positions=variant_positions,
            gene_positions=gene_positions.drop(columns="tss_position"),
        )
    with pytest.raises(ValueError, match=r"variant_positions is missing required column\(s\) \['position'\]"):
        run_eqtl_analysis(
            genotypes,
            expression,
            variant_positions=variant_positions.drop(columns="position"),
            gene_positions=gene_positions,
        )
    with pytest.raises(ValueError, match="variant_id more than once"):
        run_eqtl_analysis(
            genotypes,
            expression,
            variant_positions=pd.concat([variant_positions, variant_positions.iloc[:1]]),
            gene_positions=gene_positions,
        )
    with pytest.raises(ValueError, match="gene_id more than once"):
        run_eqtl_analysis(
            genotypes,
            expression,
            variant_positions=variant_positions,
            gene_positions=pd.concat([gene_positions, gene_positions.iloc[:1]]),
        )


def test_max_transcripts_below_one_is_still_rejected(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data

    with pytest.raises(ValueError, match="max_transcripts"):
        run_eqtl_analysis(genotypes, expression, max_transcripts=0)


def test_no_pair_in_the_window_gives_an_empty_table(data: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genotypes, expression = data
    variant_positions, gene_positions = _positions()

    result = run_eqtl_analysis(
        genotypes, expression, cis_window=0, variant_positions=variant_positions, gene_positions=gene_positions
    )

    assert result.empty
