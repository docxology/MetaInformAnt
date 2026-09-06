"""Depth tests for structural-variant multi-caller merging.

Covers consensus genotype voting, deduplication type handling, empty
callsets, and VCF-backed merging through the pysam-less fallback parser
with real VCF text written to tmp_path.
"""

from __future__ import annotations

from pathlib import Path

from metainformant.structural_variants.filtering.merge import (
    _consensus_genotype,
    _parse_vcf_basic,
    deduplicate_variants,
    merge_callsets,
    survivor_merge,
)


class TestConsensusGenotype:
    def test_missing_only_is_no_call(self) -> None:
        assert _consensus_genotype(["./.", "."]) == "./."
        assert _consensus_genotype([]) == "./."

    def test_majority_and_allele_order_normalization(self) -> None:
        # "1/0" normalizes to "0/1", so both callers agree on "0/1"
        assert _consensus_genotype(["1/0", "0/1"]) == "0/1"
        assert _consensus_genotype(["0/1", "1/1", "0/1"]) == "0/1"
        assert _consensus_genotype(["1/1", "1/1", "0/1"]) == "1/1"

    def test_merge_callsets_consensus_genotype(self) -> None:
        callsets = {
            "delly": [{"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DEL", "genotype": "1/0"}],
            "manta": [{"chrom": "chr1", "start": 1010, "end": 5010, "sv_type": "DEL", "genotype": "0/1"}],
        }
        merged, stats = merge_callsets(callsets, min_overlap=0.5)
        assert stats.n_output_variants == 1
        assert merged[0].genotype == "0/1"


class TestMergeCallsetsEdgeCases:
    def test_empty_callsets(self) -> None:
        merged, stats = merge_callsets({})
        assert merged == []
        assert stats.n_input_callsets == 0
        assert stats.n_input_variants == 0
        assert stats.n_output_variants == 0

    def test_same_caller_variants_never_merge(self) -> None:
        callsets = {
            "delly": [
                {"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DEL"},
                {"chrom": "chr1", "start": 1010, "end": 5010, "sv_type": "DEL"},
            ]
        }
        merged, stats = merge_callsets(callsets, min_overlap=0.5)
        assert stats.n_output_variants == 2
        assert all(m.n_callers == 1 for m in merged)


class TestDeduplicateVariants:
    def test_type_mismatch_keeps_both(self) -> None:
        variants = [
            {"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DEL", "quality": 60.0},
            {"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DUP", "quality": 90.0},
        ]
        deduped, n_removed = deduplicate_variants(variants, type_match=True)
        assert len(deduped) == 2
        assert n_removed == 0

    def test_type_match_false_merges_across_types(self) -> None:
        variants = [
            {"chrom": "chr1", "start": 1000, "end": 5000, "sv_type": "DEL", "quality": 60.0},
            {"chrom": "chr1", "start": 1010, "end": 5010, "sv_type": "DUP", "quality": 90.0},
        ]
        deduped, n_removed = deduplicate_variants(variants, type_match=False)
        assert len(deduped) == 1
        assert n_removed == 1
        assert deduped[0]["quality"] == 90.0


VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


class TestVcfParsing:
    def test_parse_vcf_basic_fields(self, tmp_path: Path) -> None:
        vcf = tmp_path / "calls.vcf"
        vcf.write_text(
            VCF_HEADER
            + "chr1\t1001\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=5001\n"
            + "chr1\t2001\t.\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=3001\n"
            + "# a comment line is skipped\n"
            + "chr1\tonly\tthree\tfields\n"
        )
        variants = _parse_vcf_basic(str(vcf))
        assert len(variants) == 2
        first, second = variants
        assert first["chrom"] == "chr1" and first["start"] == 1000 and first["end"] == 5001
        assert first["sv_type"] == "DEL" and first["quality"] == 50.0
        # Missing QUAL becomes 0.0
        assert second["sv_type"] == "DUP" and second["quality"] == 0.0

    def test_survivor_merge_from_vcf_files(self, tmp_path: Path) -> None:
        caller_a = tmp_path / "a.vcf"
        caller_b = tmp_path / "b.vcf"
        caller_a.write_text(VCF_HEADER + "chr1\t1001\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=5001\n")
        caller_b.write_text(VCF_HEADER + "chr1\t1021\t.\tN\t<DEL>\t70\tPASS\tSVTYPE=DEL;END=5021\n")

        merged, stats = survivor_merge([str(caller_a), str(caller_b)], max_distance=500, min_callers=2)

        assert stats.n_input_variants == 2
        assert len(merged) == 1
        assert merged[0].sv_type == "DEL"
        assert merged[0].n_callers == 2
