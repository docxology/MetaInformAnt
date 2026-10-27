"""Regression tests: ``apply_qc_filters`` validates the variant and genotype counts up front.

Defect: with ``min_qual`` unset, ``apply_qc_filters`` never compared the number of variant records with the
number of genotype columns. Fewer variant records than genotype columns surfaced as an ``IndexError`` from the
last line of the function (``vcf_data["variants"][i]``), a missing ``variants`` key as a ``KeyError``, and ragged
genotype rows as an ``IndexError`` from the transposition. Only the ``min_qual`` path checked the counts.

The expected values below are worked out by hand from the allele counts, independently of the filter code.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from metainformant.gwas.analysis.quality import apply_qc_filters

# Everything but the MAF filter is switched off, so only the minor allele frequency removes a variant.
MAF_ONLY = {"min_maf": 0.1, "max_missing": 1.0, "min_hwe_p": 0.0}

# Sample-major genotypes: GENOTYPES[s][v] is the dosage of sample s at variant v.
#   rsA  alt alleles 0+1+1+0 = 2 of 8   -> MAF 0.25   kept
#   rsB  alt alleles 1+0+0+0 = 1 of 8   -> MAF 0.125  kept
#   rsC  alt alleles 0 of 8             -> MAF 0      removed (below 0.1)
#   rsD  alt alleles 2+1+0+1 = 4 of 8   -> MAF 0.5    kept
GENOTYPES = [
    [0, 1, 0, 2],
    [1, 0, 0, 1],
    [1, 0, 0, 0],
    [0, 0, 0, 1],
]
VARIANT_IDS = ["rsA", "rsB", "rsC", "rsD"]
KEPT_IDS = ["rsA", "rsB", "rsD"]
KEPT_GENOTYPES = [[0, 1, 2], [1, 0, 1], [1, 0, 0], [0, 0, 1]]


def _variants(ids: list[str], qual: float | None = 50.0) -> list[dict[str, Any]]:
    return [
        {"chrom": "chr1", "pos": 100 * (i + 1), "id": vid, "ref": "A", "alt": ["G"], "qual": qual}
        for i, vid in enumerate(ids)
    ]


def _vcf(genotypes: Any, variants: list[dict[str, Any]] | None) -> dict[str, Any]:
    data: dict[str, Any] = {"samples": [f"S{i}" for i in range(len(genotypes))], "genotypes": genotypes}
    if variants is not None:
        data["variants"] = variants
    return data


def test_consistent_counts_still_filter_correctly() -> None:
    """Validation must not reject good data: the hand-computed survivors come out, genotypes aligned."""
    result = apply_qc_filters(_vcf(GENOTYPES, _variants(VARIANT_IDS)), MAF_ONLY)

    assert result["num_variants_before"] == 4
    assert result["num_variants_after"] == 3
    assert [v["id"] for v in result["filtered_data"]["variants"]] == KEPT_IDS
    assert np.asarray(result["filtered_data"]["genotypes"]).tolist() == KEPT_GENOTYPES


def test_numpy_genotypes_are_accepted() -> None:
    """``parse_vcf_full`` returns an ndarray; the validation has to cope with it."""
    result = apply_qc_filters(_vcf(np.array(GENOTYPES), _variants(VARIANT_IDS)), MAF_ONLY)

    assert [v["id"] for v in result["filtered_data"]["variants"]] == KEPT_IDS


def test_fewer_variant_records_than_genotype_columns_is_a_value_error() -> None:
    """Was an IndexError at the end of the function (min_qual unset)."""
    data = _vcf(GENOTYPES, _variants(VARIANT_IDS[:3]))

    with pytest.raises(ValueError, match=r"3 variant records.*4 genotype columns"):
        apply_qc_filters(data, MAF_ONLY)


def test_more_variant_records_than_genotype_columns_is_a_value_error() -> None:
    """Was accepted silently: the surplus record was dropped without a word."""
    data = _vcf(GENOTYPES, _variants(VARIANT_IDS + ["rsE"]))

    with pytest.raises(ValueError, match=r"5 variant records.*4 genotype columns"):
        apply_qc_filters(data, MAF_ONLY)


def test_missing_variants_key_is_a_value_error() -> None:
    """Was a KeyError."""
    with pytest.raises(ValueError, match=r"0 variant records.*4 genotype columns"):
        apply_qc_filters(_vcf(GENOTYPES, None), MAF_ONLY)


def test_ragged_genotype_rows_are_a_value_error() -> None:
    """A sample with fewer columns than the others was an IndexError from the transposition."""
    ragged = [row[:] for row in GENOTYPES]
    ragged[2] = ragged[2][:3]

    with pytest.raises(ValueError, match="same number of variants"):
        apply_qc_filters(_vcf(ragged, _variants(VARIANT_IDS)), MAF_ONLY)


def test_the_message_is_the_same_with_and_without_min_qual() -> None:
    """The count check no longer depends on the QUAL filter being enabled."""
    data = _vcf(GENOTYPES, _variants(VARIANT_IDS[:3]))

    with pytest.raises(ValueError, match=r"3 variant records.*4 genotype columns"):
        apply_qc_filters(data, {**MAF_ONLY, "min_qual": 30.0})


def test_data_without_genotypes_is_still_a_clean_no_op() -> None:
    """No genotype matrix at all (a sites-only VCF) keeps reporting success with nothing to filter."""
    result = apply_qc_filters({"samples": [], "genotypes": [], "variants": _variants(VARIANT_IDS)}, MAF_ONLY)

    assert result["status"] == "success"
    assert result["num_variants_before"] == 0
    assert result["num_variants_after"] == 0
