"""Regression tests for metainformant.dna.variation fixes (real computation, no test doubles)."""

import random

import pytest

from metainformant.dna.variation import calling, mutations, variants


def make_variant(chrom="chr1", pos=100, ref="A", alt=("G",), vid="rs1", qual=50.0, samples=None):
    """Build a minimal parsed-VCF variant dict."""
    return {
        "chrom": chrom,
        "pos": pos,
        "id": vid,
        "ref": ref,
        "alt": list(alt),
        "qual": qual,
        "info": {"DP": 30},
        "samples": samples or {},
    }


class TestFilterVariantsByMaf:
    def test_all_homozygous_alt_has_maf_zero_and_is_filtered(self) -> None:
        # Two hom-alt samples: alt_af = 4/4 = 1.0, so MAF must be 0.0, not 1.0
        vcf = {"variants": [make_variant(samples={"s1": {"GT": "1/1"}, "s2": {"GT": "1/1"}})]}
        result = variants.filter_variants_by_maf(vcf, min_maf=0.01)
        assert result["variants"] == []
        assert result["total_variants"] == 0

    def test_het_only_variant_kept_with_maf_half(self) -> None:
        vcf = {"variants": [make_variant(samples={"s1": {"GT": "0/1"}})]}
        result = variants.filter_variants_by_maf(vcf, min_maf=0.5)
        assert len(result["variants"]) == 1
        # min_maf above MAF filters it out
        result = variants.filter_variants_by_maf(vcf, min_maf=0.51)
        assert result["variants"] == []

    def test_missing_genotype_excluded_from_denominator(self) -> None:
        # One het call plus one missing genotype: alt_af = 1/2 = 0.5 over called genotypes
        vcf = {"variants": [make_variant(samples={"s1": {"GT": "0/1"}, "s2": {"GT": "./."}})]}
        result = variants.filter_variants_by_maf(vcf, min_maf=0.5)
        assert len(result["variants"]) == 1
        # Under the old denominator (2 samples * 2 = 4) MAF would be 0.25 and it would drop
        result = variants.filter_variants_by_maf(vcf, min_maf=0.4)
        assert len(result["variants"]) == 1

    def test_hom_ref_variant_filtered(self) -> None:
        vcf = {"variants": [make_variant(samples={"s1": {"GT": "0/0"}, "s2": {"GT": "0/0"}})]}
        result = variants.filter_variants_by_maf(vcf, min_maf=0.01)
        assert result["variants"] == []


class TestAnnotateVariants:
    def test_input_vcf_data_not_mutated(self) -> None:
        vcf = {"variants": [make_variant(vid="rs1")], "samples": ["s1"], "metadata": {}}
        annotations = {"rs1": {"GENE": "BRCA1", "CONSEQUENCE": "missense"}}
        result = variants.annotate_variants(vcf, annotations)

        # Returned data carries the annotations
        assert result["variants"][0]["info"]["GENE"] == "BRCA1"
        assert result["variants"][0]["info"]["CONSEQUENCE"] == "missense"

        # Caller's dicts are untouched
        assert vcf["variants"][0]["info"] == {"DP": 30}
        assert vcf["variants"][0]["id"] == "rs1"

    def test_unannotated_variant_shared_but_untouched(self) -> None:
        vcf = {"variants": [make_variant(vid="rs1"), make_variant(vid="rs2")]}
        annotations = {"rs1": {"GENE": "TP53"}}
        result = variants.annotate_variants(vcf, annotations)
        assert result["variants"][1] is vcf["variants"][1]
        assert "GENE" not in result["variants"][1]["info"]
        assert vcf["variants"][0]["info"] == {"DP": 30}


class TestMnvClassification:
    def test_calculate_variant_statistics_counts_mnv(self) -> None:
        vcf = {
            "variants": [
                make_variant(pos=1, ref="A", alt=["C"]),  # SNP (transversion)
                make_variant(pos=2, ref="AT", alt=["GC"]),  # MNV (equal length, multi-base)
                make_variant(pos=3, ref="AT", alt=["A"]),  # deletion
            ]
        }
        stats = variants.calculate_variant_statistics(vcf)
        assert stats["snps"] == 1
        assert stats["mnvs"] == 1
        assert stats["indels"] == 1

    def test_transition_transversion_ratio_inf_when_only_transitions(self) -> None:
        # Only A->G transitions (no transversions): ratio must be inf, not ZeroDivisionError
        vcf = {"variants": [make_variant(pos=1, ref="A", alt=["G"])]}
        stats = variants.calculate_variant_statistics(vcf)
        assert stats["transitions"] == 1
        assert stats["transversions"] == 0
        assert stats["transition_transversion_ratio"] == float("inf")

    def test_summarize_variants_by_chromosome_counts_mnv(self) -> None:
        vcf = {
            "variants": [
                make_variant(chrom="chr1", pos=1, ref="A", alt=["G"]),
                make_variant(chrom="chr1", pos=2, ref="AT", alt=["GC"]),
                make_variant(chrom="chr1", pos=3, ref="ATG", alt=["A"]),
                make_variant(chrom="chr2", pos=1, ref="C", alt=["CA"]),
            ]
        }
        summary = variants.summarize_variants_by_chromosome(vcf)
        assert summary["chr1"] == {"total": 3, "snp": 1, "mnv": 1, "indel": 1}
        assert summary["chr2"] == {"total": 1, "snp": 0, "mnv": 0, "indel": 1}


class TestPredictVariantEffect:
    coding_seq = "ATGATCGAA"

    def test_frameshift_single_base_insertion(self) -> None:
        result = variants.predict_variant_effect("A", "GG", 3, self.coding_seq)
        assert result["effect_type"] == "frameshift"

    def test_inframe_deletion(self) -> None:
        # Deleting 3 bases (ATG -> A) keeps the reading frame
        result = variants.predict_variant_effect("ATGA", "A", 3, self.coding_seq)
        assert result["effect_type"] == "inframe_indel"

    def test_inframe_insertion(self) -> None:
        # Inserting 3 bases keeps the reading frame
        result = variants.predict_variant_effect("A", "ATGA", 3, self.coding_seq)
        assert result["effect_type"] == "inframe_indel"

    def test_equal_length_multi_base_substitution_is_mnv(self) -> None:
        result = variants.predict_variant_effect("AT", "GC", 3, self.coding_seq)
        assert result["effect_type"] == "mnv"

    def test_snp_effect_unchanged(self) -> None:
        result = variants.predict_variant_effect("A", "G", 0, self.coding_seq)
        assert result["effect_type"] == "missense"
        assert result["original_aa"] == "M"
        assert result["mutated_aa"] == "V"


class TestAnnotateVariantContext:
    reference = "ACGTACGTAC"  # idx: 0 A, 1 C, 2 G, 3 T, 4 A, 5 C, 6 G, 7 T, 8 A, 9 C

    def test_snp_context_uses_one_based_position(self) -> None:
        variant = {"chrom": "chr1", "pos": 5, "ref": "A", "alt": "G"}
        annotated = calling.annotate_variant_context([variant], self.reference, window=2)
        result = annotated[0]
        # pos 5 is 1-based -> 0-based index 4, reference base 'A'
        assert result["trinucleotide_context"] == "TAC"
        assert result["trinucleotide_context"][1] == self.reference[4]
        assert result["upstream_context"] == "GT"
        assert result["downstream_context"] == "CG"
        assert result["full_context"] == "GT[A/G]CG"

    def test_sbs_channel_pyrimidine_context(self) -> None:
        # pos 2 (1-based) -> 0-based index 1, reference base 'C'; trinuc ACG
        variant = {"chrom": "chr1", "pos": 2, "ref": "C", "alt": "T"}
        result = calling.annotate_variant_context([variant], self.reference, window=1)[0]
        assert result["trinucleotide_context"] == "ACG"
        assert result["sbs_channel"] == "A[C>T]G"

    def test_context_matches_call_variants_pileup_output(self) -> None:
        pileup = [
            {
                "chrom": "chr1",
                "pos": 5,  # 1-based, as produced by call_variants_pileup
                "ref": "A",
                "depth": 30,
                "bases": {"A": 10, "G": 20},
                "quals": {"G": 35.0},
            }
        ]
        called = calling.call_variants_pileup(pileup, min_depth=10, min_qual=1.0, min_alt_freq=0.2)
        assert len(called) == 1
        annotated = calling.annotate_variant_context(called, self.reference, window=2)
        result = annotated[0]
        assert result["pos"] == 5
        assert result["trinucleotide_context"][1] == self.reference[result["pos"] - 1] == "A"


class TestGeneratePointMutations:
    def test_exact_number_of_distinct_positions_mutated(self) -> None:
        random.seed(42)
        seq = "A" * 50
        mutated = mutations.generate_point_mutations(seq, num_mutations=10)
        assert len(mutated) == len(seq)
        diffs = [i for i, (a, b) in enumerate(zip(seq, mutated)) if a != b]
        assert len(diffs) == 10
        assert all(mutated[i] != "A" for i in diffs)
        assert all(mutated[i] in "CGT" for i in diffs)

    def test_deterministic_under_seed(self) -> None:
        random.seed(7)
        first = mutations.generate_point_mutations("ACGTACGTACGT", num_mutations=4)
        random.seed(7)
        second = mutations.generate_point_mutations("ACGTACGTACGT", num_mutations=4)
        assert first == second

    def test_zero_mutations_returns_sequence(self) -> None:
        random.seed(1)
        assert mutations.generate_point_mutations("ACGTACGT", num_mutations=0) == "ACGTACGT"

    def test_too_many_mutations_raises(self) -> None:
        with pytest.raises(ValueError):
            mutations.generate_point_mutations("ACGT", num_mutations=5)


class TestCalculateSubstitutionMatrix:
    def test_counts_differences(self) -> None:
        matrix = mutations.calculate_substitution_matrix("ATCG", "AGCT")
        assert matrix == {("T", "G"): 1, ("G", "T"): 1}

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="equal length"):
            mutations.calculate_substitution_matrix("ATCG", "ATCGAA")
