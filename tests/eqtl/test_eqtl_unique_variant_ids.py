"""Regression tests: synthetic variant ids are unique, and a cis scan tests each pair once on its own genotypes.

Defect: ``create_synthetic_genotypes`` built variant ids from ``gene_id[:20]``. Kallisto transcript ids such as
``lcl|NC_037638.1_mrna_XM_623972.6_1`` share their first 20 characters (``lcl|NC_037638.1_mrna``) across a whole
chromosome, so two transcripts on one chromosome produced the same variant ids. ``cis_eqtl_scan`` then found each
duplicated id twice in a gene's cis window (testing the pair twice, under two distances) and, because
``geno.loc[id]`` returned two rows, silently took the first one: every later gene on the chromosome was tested
against the first gene's genotype vector.

Oracles: the p-values and effects of every expected test are computed with ``scipy.stats.linregress`` from the
genotype row at the same *position* in the generated matrix, which does not depend on the ids at all.
"""

from __future__ import annotations

from typing import List

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from metainformant.eqtl.workflow.synthetic import create_synthetic_genotypes, parse_gene_positions
from metainformant.gwas.finemapping.eqtl import cis_eqtl_scan

# Transcripts as Amalgkit/kallisto write them: the first two are on one chromosome, the third on another.
IDS = [
    "lcl|NC_037638.1_mrna_XM_623972.6_1",
    "lcl|NC_037638.1_mrna_XM_393171.7_2",
    "lcl|NC_037639.1_mrna_XM_624220.5_3",
]
SAMPLES = [f"SRR{9000000 + i}" for i in range(16)]
WINDOW = 500_000


def expression() -> pd.DataFrame:
    rng = np.random.default_rng(11)
    return pd.DataFrame(rng.normal(10.0, 2.0, size=(len(IDS), len(SAMPLES))), index=IDS, columns=SAMPLES)


def generated() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gene_positions = parse_gene_positions(IDS)
    genotypes, variant_positions = create_synthetic_genotypes(SAMPLES, gene_positions, variants_per_gene=3)
    return gene_positions, genotypes, variant_positions


# ---------------------------------------------------------------------------
# The generator
# ---------------------------------------------------------------------------


def test_variant_ids_of_transcripts_on_one_chromosome_are_distinct() -> None:
    _, genotypes, variant_positions = generated()

    assert len(variant_positions) == 9
    assert variant_positions["variant_id"].is_unique
    assert genotypes.index.is_unique
    assert list(genotypes.index) == list(variant_positions["variant_id"])


def test_every_variant_id_carries_the_full_gene_id() -> None:
    _, _, variant_positions = generated()

    expected = {f"var_{gene}_{i}" for gene in IDS for i in range(3)}
    assert set(variant_positions["variant_id"]) == expected


@pytest.mark.parametrize("gene_id", ["LOC1", "LOC100123", "ABCDEFGHIJKLMNOPQRST"])  # up to exactly 20 characters
def test_short_gene_ids_keep_the_variant_ids_they_always_had(gene_id: str) -> None:
    gene_positions = pd.DataFrame({"gene_id": [gene_id], "chrom": ["1"], "tss_position": [1_000_000]})

    _, variant_positions = create_synthetic_genotypes(SAMPLES, gene_positions, variants_per_gene=2)

    assert list(variant_positions["variant_id"]) == [f"var_{gene_id}_0", f"var_{gene_id}_1"]


def test_a_repeated_gene_id_is_rejected_because_it_would_repeat_variant_ids() -> None:
    gene_positions = pd.DataFrame(
        {"gene_id": ["LOC1", "LOC2", "LOC1"], "chrom": ["1", "1", "2"], "tss_position": [1, 2, 3]}
    )

    with pytest.raises(ValueError, match=r"repeats 1 gene id.*'LOC1'"):
        create_synthetic_genotypes(SAMPLES, gene_positions)


# ---------------------------------------------------------------------------
# The scan on realistic ids
# ---------------------------------------------------------------------------


def expected_tests() -> List[tuple[str, float, float]]:
    """(gene, beta, p-value) of every cis test: each gene against each variant on its chromosome in the window."""
    expr = expression()
    gene_positions, genotypes, variant_positions = generated()
    expected = []
    for _, gene in gene_positions.iterrows():
        for k, variant in variant_positions.reset_index(drop=True).iterrows():
            same_chromosome = str(variant["chrom"]) == str(gene["chrom"])
            if same_chromosome and abs(variant["position"] - gene["tss_position"]) <= WINDOW:
                fit = stats.linregress(genotypes.iloc[k].to_numpy(dtype=float), expr.loc[gene["gene_id"]].to_numpy())
                expected.append((gene["gene_id"], float(fit.slope), float(fit.pvalue)))
    return expected


def test_the_scan_tests_each_gene_variant_pair_once() -> None:
    gene_positions, genotypes, variant_positions = generated()

    scan = cis_eqtl_scan(expression(), genotypes, gene_positions, variant_positions, WINDOW, maf_threshold=0.0)

    assert not scan.duplicated(["gene_id", "variant_id"]).any()
    # genes 1 and 2 share a chromosome and each sees all six variants there; gene 3 sees its own three
    assert scan["gene_id"].value_counts().to_dict() == {IDS[0]: 6, IDS[1]: 6, IDS[2]: 3}


def test_every_scan_result_uses_the_genotypes_of_its_own_variant() -> None:
    gene_positions, genotypes, variant_positions = generated()

    scan = cis_eqtl_scan(expression(), genotypes, gene_positions, variant_positions, WINDOW, maf_threshold=0.0)

    expected = expected_tests()
    assert len(scan) == len(expected) == 15
    for gene in IDS:
        got = scan[scan["gene_id"] == gene].sort_values("pvalue")
        want = sorted((p, b) for g, b, p in expected if g == gene)
        assert got["pvalue"].tolist() == pytest.approx([p for p, _ in want], rel=1e-9)
        assert got["beta"].tolist() == pytest.approx([b for _, b in want], rel=1e-9)


def test_a_gene_is_tested_against_the_variants_of_its_neighbour_on_the_same_chromosome() -> None:
    """The cis window is positional: gene 1 reaches the three variants of gene 2, 10 kb downstream of it."""
    gene_positions, genotypes, variant_positions = generated()

    scan = cis_eqtl_scan(expression(), genotypes, gene_positions, variant_positions, WINDOW, maf_threshold=0.0)

    neighbour_variants = {f"var_{IDS[1]}_{i}" for i in range(3)}
    tested_by_gene_1 = set(scan.loc[scan["gene_id"] == IDS[0], "variant_id"])
    assert neighbour_variants <= tested_by_gene_1
    assert not any(v.startswith(f"var_{IDS[2]}") for v in tested_by_gene_1)  # another chromosome


# ---------------------------------------------------------------------------
# A repeated id in the inputs is an error, not a silent reuse
# ---------------------------------------------------------------------------


def small_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    samples = [f"s{i}" for i in range(12)]
    rng = np.random.default_rng(3)
    expr = pd.DataFrame(rng.normal(size=(2, 12)), index=["g1", "g2"], columns=samples)
    geno = pd.DataFrame(rng.integers(0, 3, size=(3, 12)), index=["v1", "v2", "v3"], columns=samples)
    gene_pos = pd.DataFrame({"gene_id": ["g1", "g2"], "chrom": ["1", "1"], "tss_position": [100_000, 200_000]})
    var_pos = pd.DataFrame({"variant_id": ["v1", "v2", "v3"], "chrom": ["1", "1", "1"], "position": [1, 2, 3]})
    return expr, geno, gene_pos, var_pos


def test_the_scan_runs_on_unique_ids() -> None:
    expr, geno, gene_pos, var_pos = small_inputs()

    scan = cis_eqtl_scan(expr, geno, gene_pos, var_pos, maf_threshold=0.0)

    assert len(scan) == 6


def test_a_repeated_variant_id_in_the_genotype_matrix_is_rejected() -> None:
    expr, geno, gene_pos, var_pos = small_inputs()
    geno.index = ["v1", "v2", "v1"]

    with pytest.raises(ValueError, match=r"genotype_matrix index has 1 repeated id.*'v1'"):
        cis_eqtl_scan(expr, geno, gene_pos, var_pos, maf_threshold=0.0)


def test_a_repeated_variant_id_in_the_variant_positions_is_rejected() -> None:
    expr, geno, gene_pos, var_pos = small_inputs()
    var_pos.loc[2, "variant_id"] = "v1"

    with pytest.raises(ValueError, match=r"variant_positions has 1 repeated id.*'v1'"):
        cis_eqtl_scan(expr, geno, gene_pos, var_pos, maf_threshold=0.0)


def test_a_repeated_gene_id_in_the_gene_positions_is_rejected() -> None:
    expr, geno, gene_pos, var_pos = small_inputs()
    gene_pos.loc[1, "gene_id"] = "g1"

    with pytest.raises(ValueError, match=r"gene_positions has 1 repeated id.*'g1'"):
        cis_eqtl_scan(expr, geno, gene_pos, var_pos, maf_threshold=0.0)


def test_a_repeated_id_is_rejected_even_when_there_are_too_few_samples_to_scan() -> None:
    """Validation comes first: bad ids are reported whatever the sample count."""
    expr, geno, gene_pos, var_pos = small_inputs()
    geno.index = ["v1", "v1", "v3"]

    with pytest.raises(ValueError, match="repeated id"):
        cis_eqtl_scan(expr.iloc[:, :4], geno.iloc[:, :4], gene_pos, var_pos, maf_threshold=0.0)
