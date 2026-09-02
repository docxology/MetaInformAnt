"""Depth tests for structural-variant quality filtering.

Covers environment-override behavior, evidence-object support extraction,
size/frequency/blacklist filter edge cases, and the population-lookup
matching internals with real variant dictionaries.
"""

from __future__ import annotations

import pytest

from metainformant.structural_variants.filtering.quality_filter import (
    _build_population_lookup,
    _lookup_population_frequency,
    apply_blacklist,
    filter_by_frequency,
    filter_by_quality,
    filter_by_size,
)


def _variant(chrom="chr1", start=1000, end=3000, sv_type="DEL", qual=50.0, support=10, **extra):
    v = {"chrom": chrom, "start": start, "end": end, "sv_type": sv_type, "quality": qual, "support": support}
    v.update(extra)
    return v


class TestFilterByQuality:
    def test_qual_and_support_thresholds(self) -> None:
        variants = [
            _variant(qual=50.0, support=10),
            _variant(qual=10.0, support=10),  # qual too low
            _variant(qual=50.0, support=1),  # too little support
        ]
        passed, stats = filter_by_quality(variants, min_qual=20.0, min_support=3)
        assert len(passed) == 1 and stats.output_count == 1
        assert stats.pass_rate == pytest.approx(1 / 3)
        assert stats.parameters["min_qual"] == 20.0

    def test_env_overrides(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("SV_MIN_QUAL", "100")
        passed, _ = filter_by_quality([_variant(qual=50.0), _variant(qual=150.0)])
        assert [v["quality"] for v in passed] == [150.0]
        monkeypatch.setenv("SV_MIN_SUPPORT", "20")
        passed, _ = filter_by_quality([_variant(qual=150.0, support=5), _variant(qual=150.0, support=25)])
        assert len(passed) == 1

    def test_support_from_evidence_object(self) -> None:
        class Evidence:
            split_reads = 6
            discordant_pairs = 4

        passed, _ = filter_by_quality([_variant(support=0, evidence=Evidence())], min_support=5)
        assert len(passed) == 1

    def test_support_from_evidence_dict(self) -> None:
        passed, _ = filter_by_quality(
            [_variant(support=0, evidence={"split_reads": 2, "discordant_pairs": 3})], min_support=4
        )
        assert len(passed) == 1

    def test_mapq_and_gq_gates(self) -> None:
        variants = [_variant(mapq=5.0), _variant(gq=10.0), _variant()]
        passed, _ = filter_by_quality(variants, min_mapq=20.0, min_gq=20.0)
        assert len(passed) == 1

    def test_alternate_qual_key(self) -> None:
        passed, _ = filter_by_quality([{"qual": 30.0, "support": 5}], min_qual=20.0)
        assert len(passed) == 1

    def test_empty_input_zero_pass_rate(self) -> None:
        passed, stats = filter_by_quality([])
        assert passed == [] and stats.pass_rate == 0.0 and stats.input_count == 0


class TestFilterBySize:
    def test_min_and_max_bounds(self) -> None:
        variants = [_variant(start=1000, end=1020), _variant(start=1000, end=3000), _variant(start=1, end=200001)]
        passed, stats = filter_by_size(variants, min_size=50, max_size=100000)
        assert stats.input_count == 3 and len(passed) == 1
        assert passed[0]["start"] == 1000 and passed[0]["end"] == 3000

    def test_translocations_bypass_size(self) -> None:
        tra = _variant(sv_type="TRA", start=0, end=0)
        passed, stats = filter_by_size([tra], min_size=50)
        assert passed == [tra] and stats.filter_name == "size_filter"

    def test_explicit_size_key_preferred(self) -> None:
        v = _variant(start=0, end=10)
        v["size"] = 5000
        passed, _ = filter_by_size([v], min_size=50)
        assert len(passed) == 1


class TestFilterByFrequency:
    def test_precomputed_af_filtered_and_annotated(self) -> None:
        variants = [_variant(af=0.05), _variant(af=0.001)]
        passed, stats = filter_by_frequency(variants, max_af=0.01)
        assert [v["af"] for v in passed] == [0.001]
        assert stats.filter_name == "frequency_filter"

    def test_allele_frequency_alias(self) -> None:
        variants = [_variant(allele_frequency=0.5), _variant(allele_frequency=0.0)]
        passed, _ = filter_by_frequency(variants, max_af=0.01)
        assert len(passed) == 1

    def test_lookup_database_list_format(self) -> None:
        db = [
            {"chrom": "chr1", "start": 1000, "end": 3000, "af": 0.9, "sv_type": "DEL"},
            {"chrom": "chr2", "start": 100, "end": 200, "af": 0.9, "sv_type": "DUP"},
        ]
        variants = [_variant(), _variant(chrom="chr3", start=5, end=50)]
        passed, stats = filter_by_frequency(variants, population_db=db, max_af=0.01)
        # Common chr1 DEL is removed; novel chr3 variant passes.
        assert [v["chrom"] for v in passed] == ["chr3"]
        assert stats.filtered_count == 1

    def test_lookup_database_dict_format(self) -> None:
        db = {"chr1:1000-3000-DEL": 0.9}
        passed, _ = filter_by_frequency([_variant()], population_db=db, max_af=0.01)
        assert passed == []

    def test_no_db_and_no_af_passes_all(self) -> None:
        variants = [_variant(), _variant(chrom="chr7")]
        passed, stats = filter_by_frequency(variants, population_db=None, max_af=0.01)
        assert passed == variants and stats.filtered_count == 0
        assert stats.parameters["population_db"] == "none"

    def test_rare_population_variant_kept(self) -> None:
        db = [{"chrom": "chr1", "start": 1000, "end": 3000, "af": 0.005, "sv_type": "DEL"}]
        passed, _ = filter_by_frequency([_variant()], population_db=db, max_af=0.01)
        assert len(passed) == 1 and passed[0]["population_af"] == 0.005


class TestBuildPopulationLookupAndMatching:
    def test_lookup_list_and_dict_agree(self) -> None:
        entries = [(1000, 3000, 0.5, "DEL")]
        from_list = _build_population_lookup(
            [{"chrom": "chr1", "start": 1000, "end": 3000, "af": 0.5, "sv_type": "DEL"}]
        )
        from_dict = _build_population_lookup({"chr1:1000-3000-DEL": 0.5})
        assert from_list == {"chr1": entries}
        assert from_dict == {"chr1": entries}

    def test_dict_keys_malformed_are_skipped(self) -> None:
        lookup = _build_population_lookup({"chrX:notanumber-3000-DEL": 0.5, "junk": 0.5})
        assert lookup == {}

    def test_lookup_reciprocal_overlap_and_window(self) -> None:
        lookup = _build_population_lookup([{"chrom": "chr1", "start": 1000, "end": 3000, "af": 0.9, "sv_type": "DEL"}])
        # Exact match
        assert _lookup_population_frequency(lookup, "chr1", 1000, 3000, "DEL", 200, 0.5) == 0.9
        # Wrong SV type does not match
        assert _lookup_population_frequency(lookup, "chr1", 1000, 3000, "DUP", 200, 0.5) is None
        # Outside window does not match
        assert _lookup_population_frequency(lookup, "chr1", 5000, 7000, "DEL", 200, 0.5) is None
        # Unknown chromosome
        assert _lookup_population_frequency(lookup, "chrQ", 1000, 3000, "DEL", 200, 0.5) is None
        # Partial overlap below reciprocal threshold
        assert _lookup_population_frequency(lookup, "chr1", 2500, 6000, "DEL", 200, 0.5) is None

    def test_best_highest_af_returned(self) -> None:
        lookup = _build_population_lookup(
            [
                {"chrom": "chr1", "start": 1000, "end": 3000, "af": 0.2, "sv_type": "DEL"},
                {"chrom": "chr1", "start": 1100, "end": 2900, "af": 0.8, "sv_type": "DEL"},
            ]
        )
        assert _lookup_population_frequency(lookup, "chr1", 1000, 3000, "DEL", 200, 0.5) == 0.8


class TestApplyBlacklist:
    def test_tuple_and_dict_regions(self) -> None:
        variants = [
            _variant(start=1000, end=3000),  # fully inside blacklist
            _variant(start=100000, end=102000),
        ]
        passed, stats = apply_blacklist(variants, [("chr1", 900, 3100)])
        assert [v["start"] for v in passed] == [100000]
        assert stats.parameters["n_blacklist_regions"] == 1

        passed_d, _ = apply_blacklist(variants, [{"chrom": "chr1", "start": 900, "end": 3100}])
        assert [v["start"] for v in passed_d] == [100000]

    def test_min_overlap_fraction_boundary(self) -> None:
        # 10 bp of a 2000 bp variant overlaps: fraction 0.005
        v = _variant(start=1000, end=3000)
        kept, _ = apply_blacklist([v], [("chr1", 2990, 3000)], min_overlap_fraction=0.01)
        assert kept == [v]
        removed, _ = apply_blacklist([v], [("chr1", 2990, 3000)], min_overlap_fraction=0.0)
        assert removed == []
        # Pure adjacency (zero-bp overlap) is not an overlap
        adjacent, _ = apply_blacklist([v], [("chr1", 3000, 3010)], min_overlap_fraction=0.0)
        assert adjacent == [v]

    def test_empty_blacklist_passes_all(self) -> None:
        variants = [_variant()]
        passed, stats = apply_blacklist(variants, [])
        assert passed == variants
        assert stats.parameters["n_blacklist_regions"] == 0

    def test_other_chromosomes_untouched(self) -> None:
        v = _variant(chrom="chr9")
        passed, _ = apply_blacklist([v], [("chr1", 0, 250000000)])
        assert passed == [v]
