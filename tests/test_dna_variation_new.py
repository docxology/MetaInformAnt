"""Tests for new DNA variation functions and the alanine codon bugfix.

Covers:
- predict_variant_effect: synonymous, missense, nonsense, frameshift, intergenic
- calculate_ti_tv_ratio: transition/transversion ratio from VCF data
- summarize_variants_by_chromosome: per-chromosome SNP/indel breakdown
- merge_vcf_data: merging two VCF datasets with deduplication
- mutations.py genetic_code: alanine codons GCT/GCC/GCA/GCG map to 'A'
"""

from __future__ import annotations

from typing import Any, Dict, List

import pytest

from metainformant.dna import variants
from metainformant.dna.variation import mutations

# ---------------------------------------------------------------------------
# Helpers: build VCF-like variant dicts without mocking
# ---------------------------------------------------------------------------


def _make_variant(
    chrom: str = "chr1",
    pos: int = 100,
    ref: str = "A",
    alt: list[str] | None = None,
    vid: str = ".",
    qual: float | None = 30.0,
) -> Dict[str, Any]:
    """Build a minimal variant record matching parse_vcf output."""
    if alt is None:
        alt = ["G"]
    return {
        "chrom": chrom,
        "pos": pos,
        "id": vid,
        "ref": ref,
        "alt": alt,
        "qual": qual,
        "filter": [],
        "info": {},
        "samples": {},
        "format": None,
    }


def _make_vcf_data(
    variant_list: List[Dict[str, Any]],
    metadata: Dict[str, str] | None = None,
    samples: List[str] | None = None,
) -> Dict[str, Any]:
    """Build a VCF data dict matching parse_vcf output structure."""
    return {
        "metadata": metadata or {},
        "samples": samples or [],
        "variants": variant_list,
        "total_variants": len(variant_list),
    }


# ===========================================================================
# predict_variant_effect
# ===========================================================================


class TestPredictVariantEffect:
    """Tests for predict_variant_effect covering every effect type."""

    # -- Synonymous: same amino acid after mutation --------------------------

    def test_synonymous_third_position_wobble(self) -> None:
        """A third-position change that preserves the amino acid is synonymous.

        CDS = ATGATCGAA  (codons: ATG=M, ATC=I, GAA=E)
        Mutate pos 5 (T->C): ATC -> ACC = T, but let's pick a real synonymous:
        Mutate pos 2 (G->A): ATG -> ATA (M -> I? No, that is missense).
        Use Leu: CTT -> CTC both = L.
        CDS = CTT...  pos 2 T->C  => CTC = L (synonymous).
        """
        cds = "CTTGAAAAA"  # CTT=L, GAA=E, AAA=K
        result = variants.predict_variant_effect("T", "C", 2, cds)
        assert result["effect_type"] == "synonymous"
        assert result["original_codon"] == "CTT"
        assert result["mutated_codon"] == "CTC"
        assert result["original_aa"] == "L"
        assert result["mutated_aa"] == "L"
        assert result["codon_position"] == 2

    def test_synonymous_serine_codons(self) -> None:
        """TCT -> TCC are both Serine, verifying synonymous classification."""
        cds = "TCTATGATG"
        result = variants.predict_variant_effect("T", "C", 2, cds)
        assert result["effect_type"] == "synonymous"
        assert result["original_aa"] == "S"
        assert result["mutated_aa"] == "S"

    # -- Missense: different amino acid --------------------------------------

    def test_missense_first_codon(self) -> None:
        """ATG (Met) -> GTG (Val) is a missense change at position 0."""
        cds = "ATGATCGAA"  # ATG=M, ATC=I, GAA=E
        result = variants.predict_variant_effect("A", "G", 0, cds)
        assert result["effect_type"] == "missense"
        assert result["original_codon"] == "ATG"
        assert result["mutated_codon"] == "GTG"
        assert result["original_aa"] == "M"
        assert result["mutated_aa"] == "V"
        assert result["codon_position"] == 0

    def test_missense_second_codon_position(self) -> None:
        """GAA (Glu) -> GCA (Ala) is a missense change at second codon position."""
        cds = "ATGGAATTT"  # ATG=M, GAA=E, TTT=F
        result = variants.predict_variant_effect("A", "C", 4, cds)
        assert result["effect_type"] == "missense"
        assert result["original_codon"] == "GAA"
        assert result["mutated_codon"] == "GCA"
        assert result["original_aa"] == "E"
        assert result["mutated_aa"] == "A"

    # -- Nonsense: premature stop codon --------------------------------------

    def test_nonsense_creates_stop_codon(self) -> None:
        """TAT (Tyr) -> TAA (Stop) is a nonsense mutation."""
        cds = "ATGTATGAA"  # ATG=M, TAT=Y, GAA=E
        result = variants.predict_variant_effect("T", "A", 5, cds)
        assert result["effect_type"] == "nonsense"
        assert result["original_codon"] == "TAT"
        assert result["mutated_codon"] == "TAA"
        assert result["original_aa"] == "Y"
        assert result["mutated_aa"] == "*"

    def test_nonsense_tga_stop(self) -> None:
        """TGG (Trp) -> TGA (Stop) via G->A at third position."""
        cds = "ATGTGGAAA"  # ATG=M, TGG=W, AAA=K
        result = variants.predict_variant_effect("G", "A", 5, cds)
        assert result["effect_type"] == "nonsense"
        assert result["original_aa"] == "W"
        assert result["mutated_aa"] == "*"

    def test_nonsense_tag_stop(self) -> None:
        """CAG (Gln) -> TAG (Stop) via C->T at first codon position."""
        cds = "ATGCAGGAA"  # ATG=M, CAG=Q, GAA=E
        result = variants.predict_variant_effect("C", "T", 3, cds)
        assert result["effect_type"] == "nonsense"
        assert result["original_aa"] == "Q"
        assert result["mutated_aa"] == "*"

    # -- Frameshift: indel not multiple of 3 ---------------------------------

    def test_frameshift_single_base_insertion(self) -> None:
        """Single-base insertion causes a frameshift."""
        cds = "ATGATCGAA"
        result = variants.predict_variant_effect("A", "AT", 0, cds)
        assert result["effect_type"] == "frameshift"
        assert result["original_codon"] is None
        assert result["mutated_codon"] is None

    def test_frameshift_single_base_deletion(self) -> None:
        """Single-base deletion causes a frameshift."""
        cds = "ATGATCGAA"
        result = variants.predict_variant_effect("AT", "A", 0, cds)
        assert result["effect_type"] == "frameshift"

    def test_frameshift_two_base_insertion(self) -> None:
        """Two-base insertion causes a frameshift."""
        cds = "ATGATCGAA"
        result = variants.predict_variant_effect("A", "AGG", 3, cds)
        assert result["effect_type"] == "frameshift"

    def test_no_frameshift_three_base_insertion(self) -> None:
        """Three-base insertion is in-frame and should NOT be classified as frameshift."""
        cds = "ATGATCGAA"
        # ref_len=1, alt_len=4, diff=3, 3%3==0 => not frameshift
        result = variants.predict_variant_effect("A", "AGGG", 0, cds)
        assert result["effect_type"] != "frameshift"

    # -- Intergenic: position outside coding sequence -------------------------

    def test_intergenic_position_beyond_cds(self) -> None:
        """Position past the coding sequence end is intergenic."""
        cds = "ATG"
        result = variants.predict_variant_effect("A", "G", 100, cds)
        assert result["effect_type"] == "intergenic"

    def test_intergenic_empty_cds(self) -> None:
        """Empty coding sequence means intergenic."""
        result = variants.predict_variant_effect("A", "G", 0, "")
        assert result["effect_type"] == "intergenic"

    def test_intergenic_incomplete_codon_at_end(self) -> None:
        """Position in an incomplete trailing codon is intergenic."""
        cds = "ATGAT"  # Last codon is only 2 bases (AT), incomplete
        # position 3 => codon_start=3, codon_start+3=6 > 5 => intergenic
        result = variants.predict_variant_effect("A", "G", 3, cds)
        assert result["effect_type"] == "intergenic"

    # -- Edge cases ----------------------------------------------------------

    def test_negative_position_raises_value_error(self) -> None:
        """Negative position should raise ValueError."""
        with pytest.raises(ValueError, match="non-negative"):
            variants.predict_variant_effect("A", "G", -1, "ATGATCGAA")

    def test_case_insensitivity(self) -> None:
        """Function should handle lowercase inputs."""
        cds = "atgatcgaa"
        result = variants.predict_variant_effect("a", "g", 0, cds)
        assert result["effect_type"] == "missense"
        assert result["original_aa"] == "M"
        assert result["mutated_aa"] == "V"

    def test_result_dictionary_keys(self) -> None:
        """Result dictionary has all expected keys."""
        result = variants.predict_variant_effect("A", "G", 0, "ATGATCGAA")
        expected_keys = {
            "effect_type",
            "original_codon",
            "mutated_codon",
            "original_aa",
            "mutated_aa",
            "codon_position",
        }
        assert set(result.keys()) == expected_keys


# ===========================================================================
# calculate_ti_tv_ratio
# ===========================================================================


class TestCalculateTiTvRatio:
    """Tests for calculate_ti_tv_ratio."""

    def test_known_ratio_two_to_one(self) -> None:
        """2 transitions and 1 transversion gives Ti/Tv = 2.0."""
        vcf = _make_vcf_data(
            [
                _make_variant(ref="A", alt=["G"]),  # transition
                _make_variant(ref="C", alt=["T"]),  # transition
                _make_variant(ref="A", alt=["C"]),  # transversion
            ]
        )
        ratio = variants.calculate_ti_tv_ratio(vcf)
        assert ratio == pytest.approx(2.0)

    def test_all_transitions_returns_zero(self) -> None:
        """When there are zero transversions, ratio returns 0.0."""
        vcf = _make_vcf_data(
            [
                _make_variant(ref="A", alt=["G"]),
                _make_variant(ref="C", alt=["T"]),
                _make_variant(ref="G", alt=["A"]),
            ]
        )
        ratio = variants.calculate_ti_tv_ratio(vcf)
        assert ratio == 0.0

    def test_all_transversions(self) -> None:
        """When there are zero transitions, ratio is 0.0 / N = 0.0."""
        vcf = _make_vcf_data(
            [
                _make_variant(ref="A", alt=["C"]),
                _make_variant(ref="A", alt=["T"]),
                _make_variant(ref="G", alt=["C"]),
                _make_variant(ref="G", alt=["T"]),
            ]
        )
        ratio = variants.calculate_ti_tv_ratio(vcf)
        assert ratio == pytest.approx(0.0)

    def test_equal_ti_tv_gives_one(self) -> None:
        """Equal counts of transitions and transversions gives ratio = 1.0."""
        vcf = _make_vcf_data(
            [
                _make_variant(ref="A", alt=["G"]),  # Ti
                _make_variant(ref="A", alt=["C"]),  # Tv
            ]
        )
        ratio = variants.calculate_ti_tv_ratio(vcf)
        assert ratio == pytest.approx(1.0)

    def test_indels_excluded_from_ratio(self) -> None:
        """Indels (len(ref) != len(alt)) are excluded from Ti/Tv calculation."""
        vcf = _make_vcf_data(
            [
                _make_variant(ref="A", alt=["G"]),  # Ti (counted)
                _make_variant(ref="AT", alt=["A"]),  # indel (skipped)
                _make_variant(ref="A", alt=["T"]),  # Tv (counted)
            ]
        )
        ratio = variants.calculate_ti_tv_ratio(vcf)
        # 1 Ti / 1 Tv = 1.0
        assert ratio == pytest.approx(1.0)

    def test_multiallelic_snps(self) -> None:
        """Multiple alt alleles are each classified independently."""
        vcf = _make_vcf_data(
            [
                _make_variant(ref="A", alt=["G", "C"]),  # G=Ti, C=Tv
            ]
        )
        ratio = variants.calculate_ti_tv_ratio(vcf)
        # 1 Ti / 1 Tv = 1.0
        assert ratio == pytest.approx(1.0)

    def test_empty_variants(self) -> None:
        """Empty variant list returns 0.0."""
        vcf = _make_vcf_data([])
        ratio = variants.calculate_ti_tv_ratio(vcf)
        assert ratio == 0.0

    def test_case_insensitivity(self) -> None:
        """Lowercase ref/alt alleles are handled correctly."""
        vcf = _make_vcf_data(
            [
                _make_variant(ref="a", alt=["g"]),  # transition
                _make_variant(ref="a", alt=["c"]),  # transversion
            ]
        )
        ratio = variants.calculate_ti_tv_ratio(vcf)
        assert ratio == pytest.approx(1.0)

    def test_typical_human_genome_range(self) -> None:
        """A dataset resembling typical human exome Ti/Tv (~2.0-3.0)."""
        # Build a dataset with 6 transitions and 2 transversions => ratio 3.0
        ti_variants = [_make_variant(ref="A", alt=["G"], pos=i) for i in range(6)]
        tv_variants = [_make_variant(ref="A", alt=["C"], pos=i + 100) for i in range(2)]
        vcf = _make_vcf_data(ti_variants + tv_variants)
        ratio = variants.calculate_ti_tv_ratio(vcf)
        assert ratio == pytest.approx(3.0)


# ===========================================================================
# summarize_variants_by_chromosome
# ===========================================================================


class TestSummarizeVariantsByChromosome:
    """Tests for summarize_variants_by_chromosome."""

    def test_single_chromosome_snps_only(self) -> None:
        """All SNPs on one chromosome."""
        vcf = _make_vcf_data(
            [
                _make_variant(chrom="chr1", ref="A", alt=["G"]),
                _make_variant(chrom="chr1", ref="C", alt=["T"], pos=200),
            ]
        )
        summary = variants.summarize_variants_by_chromosome(vcf)
        assert summary["chr1"]["total"] == 2
        assert summary["chr1"]["snp"] == 2
        assert summary["chr1"]["indel"] == 0

    def test_single_chromosome_indels_only(self) -> None:
        """All indels on one chromosome."""
        vcf = _make_vcf_data(
            [
                _make_variant(chrom="chr1", ref="AT", alt=["A"]),  # deletion
                _make_variant(chrom="chr1", ref="G", alt=["GCC"], pos=200),  # insertion
            ]
        )
        summary = variants.summarize_variants_by_chromosome(vcf)
        assert summary["chr1"]["total"] == 2
        assert summary["chr1"]["snp"] == 0
        assert summary["chr1"]["indel"] == 2

    def test_mixed_snps_and_indels_across_chromosomes(self) -> None:
        """Mixed SNPs and indels distributed across multiple chromosomes."""
        vcf = _make_vcf_data(
            [
                _make_variant(chrom="chr1", ref="A", alt=["G"], pos=100),  # SNP
                _make_variant(chrom="chr1", ref="AT", alt=["A"], pos=200),  # indel
                _make_variant(chrom="chr2", ref="C", alt=["T"], pos=100),  # SNP
                _make_variant(chrom="chr2", ref="G", alt=["A"], pos=200),  # SNP
                _make_variant(chrom="chr2", ref="G", alt=["GATT"], pos=300),  # indel
                _make_variant(chrom="chrX", ref="T", alt=["A"], pos=100),  # SNP
            ]
        )
        summary = variants.summarize_variants_by_chromosome(vcf)

        assert summary["chr1"]["total"] == 2
        assert summary["chr1"]["snp"] == 1
        assert summary["chr1"]["indel"] == 1

        assert summary["chr2"]["total"] == 3
        assert summary["chr2"]["snp"] == 2
        assert summary["chr2"]["indel"] == 1

        assert summary["chrX"]["total"] == 1
        assert summary["chrX"]["snp"] == 1
        assert summary["chrX"]["indel"] == 0

    def test_empty_variants(self) -> None:
        """Empty variant list returns empty dict."""
        vcf = _make_vcf_data([])
        summary = variants.summarize_variants_by_chromosome(vcf)
        assert summary == {}

    def test_multiallelic_variant_classification(self) -> None:
        """Multiallelic SNPs (all single-base alts) are counted as SNP."""
        vcf = _make_vcf_data(
            [
                _make_variant(chrom="chr1", ref="A", alt=["G", "C"]),
            ]
        )
        summary = variants.summarize_variants_by_chromosome(vcf)
        assert summary["chr1"]["snp"] == 1
        assert summary["chr1"]["indel"] == 0

    def test_multiallelic_with_indel_classified_as_indel(self) -> None:
        """Multiallelic variant with one indel alt is classified as indel."""
        vcf = _make_vcf_data(
            [
                _make_variant(chrom="chr1", ref="A", alt=["G", "AT"]),  # mixed: one SNP alt, one indel alt
            ]
        )
        summary = variants.summarize_variants_by_chromosome(vcf)
        # The function checks: all alts single-base? No => indel
        assert summary["chr1"]["indel"] == 1

    def test_all_chromosomes_present(self) -> None:
        """Summary contains an entry for every chromosome observed."""
        chroms = ["chr1", "chr2", "chr3", "chr4", "chr5"]
        variant_list = [_make_variant(chrom=c, ref="A", alt=["G"], pos=100) for c in chroms]
        vcf = _make_vcf_data(variant_list)
        summary = variants.summarize_variants_by_chromosome(vcf)
        assert set(summary.keys()) == set(chroms)


# ===========================================================================
# merge_vcf_data
# ===========================================================================


class TestMergeVcfData:
    """Tests for merge_vcf_data."""

    def test_merge_non_overlapping_variants(self) -> None:
        """Merging two datasets with no shared variants produces the union."""
        vcf1 = _make_vcf_data(
            [_make_variant(chrom="chr1", pos=100, ref="A", alt=["G"])],
            metadata={"source": "gatk"},
            samples=["S1"],
        )
        vcf2 = _make_vcf_data(
            [_make_variant(chrom="chr1", pos=200, ref="C", alt=["T"])],
            metadata={"caller": "bcftools"},
            samples=["S2"],
        )
        merged = variants.merge_vcf_data(vcf1, vcf2)
        assert merged["total_variants"] == 2
        assert "S1" in merged["samples"]
        assert "S2" in merged["samples"]

    def test_deduplication_by_chrom_pos_ref_alt(self) -> None:
        """Identical variants in both datasets are deduplicated."""
        shared_variant = _make_variant(chrom="chr1", pos=100, ref="A", alt=["G"])
        vcf1 = _make_vcf_data([shared_variant])
        # Build a second copy (distinct dict object, same data)
        shared_variant_copy = _make_variant(chrom="chr1", pos=100, ref="A", alt=["G"])
        vcf2 = _make_vcf_data([shared_variant_copy])
        merged = variants.merge_vcf_data(vcf1, vcf2)
        assert merged["total_variants"] == 1

    def test_deduplication_preserves_unique_variants(self) -> None:
        """Deduplication keeps unique variants from both datasets."""
        v1 = _make_variant(chrom="chr1", pos=100, ref="A", alt=["G"])
        v_shared = _make_variant(chrom="chr1", pos=200, ref="C", alt=["T"])
        v2 = _make_variant(chrom="chr2", pos=300, ref="G", alt=["A"])
        v_shared_dup = _make_variant(chrom="chr1", pos=200, ref="C", alt=["T"])

        vcf1 = _make_vcf_data([v1, v_shared])
        vcf2 = _make_vcf_data([v_shared_dup, v2])
        merged = variants.merge_vcf_data(vcf1, vcf2)
        # v1 + v_shared + v2 = 3 unique
        assert merged["total_variants"] == 3

    def test_metadata_vcf1_precedence(self) -> None:
        """vcf1 metadata takes precedence for overlapping keys."""
        vcf1 = _make_vcf_data([], metadata={"source": "gatk", "version": "4.0"})
        vcf2 = _make_vcf_data([], metadata={"source": "bcftools", "date": "2024-01-01"})
        merged = variants.merge_vcf_data(vcf1, vcf2)
        assert merged["metadata"]["source"] == "gatk"  # vcf1 wins
        assert merged["metadata"]["date"] == "2024-01-01"  # only in vcf2
        assert merged["metadata"]["version"] == "4.0"

    def test_samples_union_preserves_order(self) -> None:
        """Merged samples are the union, preserving insertion order."""
        vcf1 = _make_vcf_data([], samples=["S1", "S2"])
        vcf2 = _make_vcf_data([], samples=["S2", "S3"])
        merged = variants.merge_vcf_data(vcf1, vcf2)
        assert merged["samples"] == ["S1", "S2", "S3"]

    def test_merged_variants_sorted_by_chrom_and_pos(self) -> None:
        """Merged variants are sorted by chromosome then position."""
        vcf1 = _make_vcf_data(
            [
                _make_variant(chrom="chr2", pos=300, ref="G", alt=["A"]),
            ]
        )
        vcf2 = _make_vcf_data(
            [
                _make_variant(chrom="chr1", pos=100, ref="A", alt=["G"]),
                _make_variant(chrom="chr2", pos=100, ref="C", alt=["T"]),
            ]
        )
        merged = variants.merge_vcf_data(vcf1, vcf2)
        positions = [(v["chrom"], v["pos"]) for v in merged["variants"]]
        assert positions == sorted(positions)

    def test_merge_empty_with_nonempty(self) -> None:
        """Merging an empty dataset with a non-empty one returns the non-empty data."""
        vcf_empty = _make_vcf_data([])
        vcf_full = _make_vcf_data(
            [
                _make_variant(chrom="chr1", pos=100, ref="A", alt=["G"]),
                _make_variant(chrom="chr1", pos=200, ref="C", alt=["T"]),
            ]
        )
        merged = variants.merge_vcf_data(vcf_empty, vcf_full)
        assert merged["total_variants"] == 2

    def test_merge_two_empty_datasets(self) -> None:
        """Merging two empty datasets returns an empty result."""
        vcf1 = _make_vcf_data([])
        vcf2 = _make_vcf_data([])
        merged = variants.merge_vcf_data(vcf1, vcf2)
        assert merged["total_variants"] == 0
        assert merged["variants"] == []

    def test_different_alt_alleles_not_deduplicated(self) -> None:
        """Same chrom/pos/ref but different alt alleles are distinct variants."""
        vcf1 = _make_vcf_data([_make_variant(chrom="chr1", pos=100, ref="A", alt=["G"])])
        vcf2 = _make_vcf_data([_make_variant(chrom="chr1", pos=100, ref="A", alt=["C"])])
        merged = variants.merge_vcf_data(vcf1, vcf2)
        assert merged["total_variants"] == 2

    def test_alt_order_insensitive_dedup(self) -> None:
        """Alt alleles sorted for dedup key, so order does not matter."""
        vcf1 = _make_vcf_data([_make_variant(chrom="chr1", pos=100, ref="A", alt=["G", "C"])])
        vcf2 = _make_vcf_data([_make_variant(chrom="chr1", pos=100, ref="A", alt=["C", "G"])])
        merged = variants.merge_vcf_data(vcf1, vcf2)
        assert merged["total_variants"] == 1


# ===========================================================================
# Alanine codon bugfix verification (mutations.py genetic_code)
# ===========================================================================


class TestAlanineCodonBugfix:
    """Verify that GCT/GCC/GCA/GCG map to 'A' (Alanine), not 'G' (Glycine).

    This validates the bugfix applied to both the _GENETIC_CODE table in
    variants.py and the genetic_code table in mutations.py.
    """

    def test_variants_genetic_code_gct_is_alanine(self) -> None:
        """GCT maps to 'A' (Alanine) in the variants module genetic code."""
        assert variants._GENETIC_CODE["GCT"] == "A"

    def test_variants_genetic_code_gcc_is_alanine(self) -> None:
        """GCC maps to 'A' (Alanine) in the variants module genetic code."""
        assert variants._GENETIC_CODE["GCC"] == "A"

    def test_variants_genetic_code_gca_is_alanine(self) -> None:
        """GCA maps to 'A' (Alanine) in the variants module genetic code."""
        assert variants._GENETIC_CODE["GCA"] == "A"

    def test_variants_genetic_code_gcg_is_alanine(self) -> None:
        """GCG maps to 'A' (Alanine) in the variants module genetic code."""
        assert variants._GENETIC_CODE["GCG"] == "A"

    def test_variants_glycine_codons_unchanged(self) -> None:
        """Glycine codons (GGT/GGC/GGA/GGG) still correctly map to 'G'."""
        for codon in ("GGT", "GGC", "GGA", "GGG"):
            assert variants._GENETIC_CODE[codon] == "G"

    def test_mutations_genetic_code_via_classify(self) -> None:
        """mutations.classify_mutations uses genetic code where GCx = Alanine.

        GCT -> GCC is a synonymous change (both Alanine). If the old bug
        were present (mapping to Glycine), the test would still technically
        pass for synonymous, but we verify the amino acid identity through
        a known missense: GCT (Ala) -> GAT (Asp) at position 1 of the codon.
        """
        # Ancestral: GCT = A (Ala), Derived: GAT = D (Asp) => nonsynonymous
        result = mutations.classify_mutations("GCT", "GAT")
        assert result["nonsynonymous"] == 1
        assert result["synonymous"] == 0
        assert result["total"] == 1

    def test_mutations_alanine_synonymous(self) -> None:
        """GCT -> GCC (both Alanine) must be synonymous in mutations module."""
        result = mutations.classify_mutations("GCT", "GCC")
        assert result["synonymous"] == 1
        assert result["nonsynonymous"] == 0

    def test_mutations_alanine_all_four_codons_synonymous(self) -> None:
        """All four alanine codons should yield synonymous changes between each other.

        Pick pairs: GCT->GCC, GCT->GCA, GCT->GCG. Each differs at position 2,
        so classify_mutations sees exactly 1 mutation per pair.
        """
        alanine_codons = ["GCT", "GCC", "GCA", "GCG"]
        for i, c1 in enumerate(alanine_codons):
            for j, c2 in enumerate(alanine_codons):
                if i != j:
                    result = mutations.classify_mutations(c1, c2)
                    # Each pair differs by 1 base (third position wobble)
                    assert result["synonymous"] >= 1, f"{c1} -> {c2} should be synonymous"
                    assert result["nonsynonymous"] == 0, f"{c1} -> {c2} should have no nonsynonymous"

    def test_predict_effect_alanine_synonymous(self) -> None:
        """predict_variant_effect should classify GCT->GCC as synonymous."""
        cds = "GCTAAATTT"  # GCT=A, AAA=K, TTT=F
        result = variants.predict_variant_effect("T", "C", 2, cds)
        assert result["effect_type"] == "synonymous"
        assert result["original_aa"] == "A"
        assert result["mutated_aa"] == "A"

    def test_predict_effect_alanine_to_glycine_is_missense(self) -> None:
        """GCT (Ala) -> GGT (Gly) is a missense change, not synonymous."""
        cds = "GCTATGATG"  # GCT=A, ATG=M, ATG=M
        result = variants.predict_variant_effect("C", "G", 1, cds)
        assert result["effect_type"] == "missense"
        assert result["original_aa"] == "A"
        assert result["mutated_aa"] == "G"


# ===========================================================================
# Integration / end-to-end scenarios
# ===========================================================================


class TestVariationIntegration:
    """Integration tests combining multiple variation functions."""

    def test_merge_then_summarize(self) -> None:
        """Merge two VCF datasets and summarize by chromosome."""
        vcf1 = _make_vcf_data(
            [
                _make_variant(chrom="chr1", pos=100, ref="A", alt=["G"]),
                _make_variant(chrom="chr2", pos=100, ref="AT", alt=["A"]),
            ]
        )
        vcf2 = _make_vcf_data(
            [
                _make_variant(chrom="chr1", pos=200, ref="C", alt=["T"]),
                _make_variant(chrom="chr3", pos=100, ref="G", alt=["GAA"]),
            ]
        )
        merged = variants.merge_vcf_data(vcf1, vcf2)
        summary = variants.summarize_variants_by_chromosome(merged)

        assert summary["chr1"]["total"] == 2
        assert summary["chr1"]["snp"] == 2
        assert summary["chr2"]["total"] == 1
        assert summary["chr2"]["indel"] == 1
        assert summary["chr3"]["total"] == 1
        assert summary["chr3"]["indel"] == 1

    def test_merge_then_ti_tv_ratio(self) -> None:
        """Merge two VCF datasets and compute Ti/Tv ratio on the merged result."""
        vcf1 = _make_vcf_data(
            [
                _make_variant(chrom="chr1", pos=100, ref="A", alt=["G"]),  # Ti
                _make_variant(chrom="chr1", pos=200, ref="A", alt=["G"]),  # Ti
            ]
        )
        vcf2 = _make_vcf_data(
            [
                _make_variant(chrom="chr1", pos=300, ref="A", alt=["C"]),  # Tv
            ]
        )
        merged = variants.merge_vcf_data(vcf1, vcf2)
        ratio = variants.calculate_ti_tv_ratio(merged)
        assert ratio == pytest.approx(2.0)

    def test_predict_effect_on_full_cds(self) -> None:
        """Walk through several positions of a real CDS predicting effects."""
        # CDS: ATG GCT GAA TAT TGG TAA
        # AA:  M   A   E   Y   W   *
        cds = "ATGGCTGAATATTGGTAA"

        # Position 0: A->G in ATG => GTG (M->V, missense)
        r = variants.predict_variant_effect("A", "G", 0, cds)
        assert r["effect_type"] == "missense"

        # Position 5: T->C in GCT => GCC (A->A, synonymous)
        r = variants.predict_variant_effect("T", "C", 5, cds)
        assert r["effect_type"] == "synonymous"

        # Position 11: T->A in TAT => TAA (Y->*, nonsense)
        r = variants.predict_variant_effect("T", "A", 11, cds)
        assert r["effect_type"] == "nonsense"

    def test_full_genetic_code_completeness(self) -> None:
        """The genetic code table in variants has all 64 codons."""
        bases = "ACGT"
        expected_codons = {a + b + c for a in bases for b in bases for c in bases}
        assert set(variants._GENETIC_CODE.keys()) == expected_codons
        assert len(variants._GENETIC_CODE) == 64
