"""Regression tests for ``survivor_merge`` genotypes and qualities.

``survivor_merge`` never looked at the genotype of the calls it merged, so every merged variant came out as
``./.``, and its quality was ``10 * n_callers * log10(n_callers + 1)``: the same number for two calls of quality
5 and two calls of quality 900. It now combines the callers' evidence by a documented rule:

* genotype: majority vote among the calls' non-missing genotypes, allele order ignored; a tie goes to the
  genotype with the larger summed call quality, then to the one that sorts first;
* quality: the mean call quality plus the agreement bonus ``10 * n_callers * log10(n_callers + 1)``.

VCF input carries its genotype in the first sample column, which the loaders used to drop.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pytest

from metainformant.structural_variants.filtering.merge import _parse_vcf_basic, survivor_merge


def _call(caller: str, start: int, genotype: str | None = None, quality: float | None = None, **extra: Any) -> dict:
    call: dict[str, Any] = {
        "chrom": "chr1",
        "start": start,
        "end": start + 4_000,
        "sv_type": "DEL",
        "_caller": caller,
    }
    if genotype is not None:
        call["genotype"] = genotype
    if quality is not None:
        call["quality"] = quality
    call.update(extra)
    return call


def _merge(calls: list[dict], min_callers: int = 2) -> Any:
    merged, _ = survivor_merge(calls, max_distance=500, min_callers=min_callers)
    assert len(merged) == 1
    return merged[0]


def _bonus(n_callers: int) -> float:
    return 10.0 * n_callers * math.log10(n_callers + 1)


class TestSurvivorGenotype:
    def test_agreeing_callers_give_their_genotype_whatever_the_allele_order(self) -> None:
        merged = _merge([_call("delly", 1_000, "0/1"), _call("manta", 1_010, "1/0")])
        assert merged.genotype == "0/1"

    def test_the_majority_wins(self) -> None:
        merged = _merge(
            [_call("delly", 1_000, "0/1", 10), _call("manta", 1_010, "0/1", 10), _call("lumpy", 1_020, "1/1", 900)]
        )
        assert merged.genotype == "0/1"

    def test_a_tie_goes_to_the_genotype_with_the_larger_summed_quality(self) -> None:
        assert _merge([_call("delly", 1_000, "1/1", 90), _call("manta", 1_010, "0/1", 30)]).genotype == "1/1"
        assert _merge([_call("delly", 1_000, "1/1", 30), _call("manta", 1_010, "0/1", 90)]).genotype == "0/1"

    def test_a_tie_in_votes_and_quality_goes_to_the_genotype_that_sorts_first(self) -> None:
        assert _merge([_call("delly", 1_000, "1/1", 50), _call("manta", 1_010, "0/1", 50)]).genotype == "0/1"
        assert _merge([_call("delly", 1_000, "0/1"), _call("manta", 1_010, "1/1")]).genotype == "0/1"

    def test_missing_genotypes_do_not_vote(self) -> None:
        assert _merge([_call("delly", 1_000, "./."), _call("manta", 1_010, "1/1")]).genotype == "1/1"
        assert _merge([_call("delly", 1_000), _call("manta", 1_010, "0/0")]).genotype == "0/0"

    def test_a_genotype_of_none_is_a_missing_genotype(self) -> None:
        with_none = dict(_call("delly", 1_000), genotype=None)  # the key is present, its value is None
        assert _merge([with_none, _call("manta", 1_010, "1/1")]).genotype == "1/1"
        merged = _merge([with_none, dict(_call("manta", 1_010), genotype=None)])
        assert merged.genotype == "./."

    def test_without_any_genotype_the_call_stays_missing(self) -> None:
        assert _merge([_call("delly", 1_000), _call("manta", 1_010)]).genotype == "./."
        assert _merge([_call("delly", 1_000, "."), _call("manta", 1_010, "./.")]).genotype == "./."


class TestSurvivorQuality:
    def test_quality_is_the_mean_call_quality_plus_the_agreement_bonus(self) -> None:
        merged = _merge([_call("delly", 1_000, quality=80.0), _call("manta", 1_010, quality=40.0)])
        assert merged.consensus_quality == pytest.approx(60.0 + _bonus(2))
        assert merged.consensus_quality == pytest.approx(69.5424, abs=1e-4)

    def test_quality_follows_the_quality_of_the_calls(self) -> None:
        low = _merge([_call("delly", 1_000, quality=5.0), _call("manta", 1_010, quality=5.0)])
        high = _merge([_call("delly", 1_000, quality=900.0), _call("manta", 1_010, quality=900.0)])
        assert high.consensus_quality - low.consensus_quality == pytest.approx(895.0)

    def test_more_agreeing_callers_raise_the_quality(self) -> None:
        two = _merge([_call("a", 1_000, quality=30.0), _call("b", 1_010, quality=30.0)])
        three = _merge(
            [_call("a", 1_000, quality=30.0), _call("b", 1_010, quality=30.0), _call("c", 1_020, quality=30.0)]
        )
        assert three.consensus_quality == pytest.approx(30.0 + _bonus(3))
        assert three.consensus_quality > two.consensus_quality

    def test_single_caller_groups_use_their_own_quality(self) -> None:
        merged = _merge([_call("delly", 1_000, quality=33.0)], min_callers=1)
        assert merged.consensus_quality == pytest.approx(33.0 + _bonus(1))

    def test_calls_without_a_quality_score_the_agreement_bonus_alone(self) -> None:
        merged = _merge([_call("delly", 1_000), _call("manta", 1_010)])
        assert merged.consensus_quality == pytest.approx(_bonus(2))

    def test_qual_is_an_alias_and_none_counts_as_zero(self) -> None:
        merged = _merge([_call("delly", 1_000, qual=70.0), _call("manta", 1_010, quality=None, qual=None)])
        assert merged.consensus_quality == pytest.approx(35.0 + _bonus(2))

    def test_the_merged_variant_keeps_its_callers_in_a_stable_order(self) -> None:
        merged = _merge([_call("manta", 1_000), _call("delly", 1_010), _call("lumpy", 1_020)])
        assert merged.caller_names == ["delly", "lumpy", "manta"]
        assert merged.n_callers == 3


# ---------------------------------------------------------------------------
# VCF input
# ---------------------------------------------------------------------------

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"


def _write_vcf(path: Path, *records: str, header: str = HEADER) -> str:
    path.write_text(header + "".join(record + "\n" for record in records))
    return str(path)


class TestVcfGenotypes:
    def test_genotype_of_the_first_sample_is_loaded(self, tmp_path: Path) -> None:
        vcf = _write_vcf(tmp_path / "a.vcf", "chr1\t1001\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=5001\tGT:GQ\t0/1:30")
        [record] = _parse_vcf_basic(vcf)
        assert record["genotype"] == "0/1"

    def test_genotype_is_found_wherever_gt_sits_in_format(self, tmp_path: Path) -> None:
        vcf = _write_vcf(tmp_path / "a.vcf", "chr1\t1001\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=5001\tGQ:GT\t30:1/1")
        assert _parse_vcf_basic(vcf)[0]["genotype"] == "1/1"

    def test_phased_and_missing_genotypes_are_written_with_slashes(self, tmp_path: Path) -> None:
        vcf = _write_vcf(
            tmp_path / "a.vcf",
            "chr1\t1001\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=5001\tGT\t0|1",
            "chr1\t2001\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=6001\tGT\t.",
            "chr1\t3001\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=7001\tGT\t.|.",
        )
        assert [r["genotype"] for r in _parse_vcf_basic(vcf)] == ["0/1", "./.", "./."]

    def test_records_without_a_genotype_have_no_genotype_key(self, tmp_path: Path) -> None:
        sites_only = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        vcf = _write_vcf(
            tmp_path / "a.vcf", "chr1\t1001\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=5001", header=sites_only
        )
        assert "genotype" not in _parse_vcf_basic(vcf)[0]
        no_gt = _write_vcf(tmp_path / "b.vcf", "chr1\t1001\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=5001\tGQ\t30")
        assert "genotype" not in _parse_vcf_basic(no_gt)[0]

    def test_survivor_merge_reads_genotypes_from_vcf_files(self, tmp_path: Path) -> None:
        first = _write_vcf(
            tmp_path / "delly.vcf", "chr1\t1001\t.\tN\t<DEL>\t80\tPASS\tSVTYPE=DEL;END=5001\tGT:GQ\t0/1:30"
        )
        second = _write_vcf(
            tmp_path / "manta.vcf", "chr1\t1021\t.\tN\t<DEL>\t40\tPASS\tSVTYPE=DEL;END=5021\tGT:GQ\t0|1:25"
        )

        merged, stats = survivor_merge([first, second], max_distance=500, min_callers=2)

        assert stats.n_input_variants == 2
        assert len(merged) == 1
        assert merged[0].genotype == "0/1"
        assert merged[0].consensus_quality == pytest.approx(60.0 + _bonus(2))

    def test_survivor_merge_breaks_a_genotype_tie_by_vcf_quality(self, tmp_path: Path) -> None:
        first = _write_vcf(tmp_path / "a.vcf", "chr1\t1001\t.\tN\t<DEL>\t10\tPASS\tSVTYPE=DEL;END=5001\tGT\t0/1")
        second = _write_vcf(tmp_path / "b.vcf", "chr1\t1021\t.\tN\t<DEL>\t90\tPASS\tSVTYPE=DEL;END=5021\tGT\t1/1")
        merged, _ = survivor_merge([first, second], max_distance=500, min_callers=2)
        assert merged[0].genotype == "1/1"


class TestGenotypeHelpers:
    @pytest.mark.parametrize(
        ("alleles", "expected"),
        [
            ((0, 1), "0/1"),
            ((1, 1), "1/1"),
            ((None, None), "./."),
            ((0, None), "0/."),
            ((1,), "1"),
            ((), None),
            (None, None),
        ],
    )
    def test_format_genotype(self, alleles: Any, expected: str | None) -> None:
        from metainformant.structural_variants.filtering.merge import _format_genotype

        assert _format_genotype(alleles) == expected

    @pytest.mark.parametrize(
        ("fields", "expected"),
        [
            (["chr1", "1", ".", "N", "<DEL>", "5", "PASS", "SVTYPE=DEL", "GT", "0/1"], "0/1"),
            (["chr1", "1", ".", "N", "<DEL>", "5", "PASS", "SVTYPE=DEL", "GT:GQ", "1|1:9"], "1/1"),
            (["chr1", "1", ".", "N", "<DEL>", "5", "PASS", "SVTYPE=DEL", "GQ:GT", "9:0/0"], "0/0"),
            (["chr1", "1", ".", "N", "<DEL>", "5", "PASS", "SVTYPE=DEL", "GT", "."], "./."),
            (["chr1", "1", ".", "N", "<DEL>", "5", "PASS", "SVTYPE=DEL", "GQ", "9"], None),
            (["chr1", "1", ".", "N", "<DEL>", "5", "PASS", "SVTYPE=DEL", "GQ:GT", "9"], None),
            (["chr1", "1", ".", "N", "<DEL>", "5", "PASS", "SVTYPE=DEL"], None),
        ],
    )
    def test_vcf_sample_genotype(self, fields: list[str], expected: str | None) -> None:
        from metainformant.structural_variants.filtering.merge import _vcf_sample_genotype

        assert _vcf_sample_genotype(fields) == expected
