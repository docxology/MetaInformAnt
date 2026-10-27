"""Tests for metagenomics shotgun submodule.

Tests assembly, binning, community profiling.
Uses real implementations -- real-implementation policy.
"""

from __future__ import annotations

import os
import random
import subprocess
import sys
import textwrap

import pytest

from metainformant.metagenomics.shotgun.assembly import (
    AssemblyStats,
    Contig,
    Scaffold,
    assemble_contigs,
    calculate_assembly_stats,
    scaffold_contigs,
)
from metainformant.metagenomics.shotgun.binning import (
    BinningResult,
    GenomeBin,
    assess_bin_quality,
    bin_contigs,
    calculate_tetranucleotide_freq,
    refine_bins,
)
from metainformant.metagenomics.shotgun.profiling import (
    CommunityProfile,
    KmerIndex,
    TaxonProfile,
    build_kmer_index,
    calculate_relative_abundance,
    profile_community,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Ten expected single-copy marker genes used by the bin-quality tests
EXPECTED_MARKERS = [f"m{i}" for i in range(1, 11)]


@pytest.fixture
def overlapping_reads() -> list[str]:
    """Reads that share a common region for assembly."""
    base = "ATCGATCGATCGATCGATCGATCGATCGATCG"
    return [
        base,
        base,
        base[5:] + "AAAAA",
        "TTTTT" + base[:27],
        base,
    ]


@pytest.fixture
def simple_contigs() -> list[Contig]:
    """Pre-made contigs for scaffold/stats tests."""
    return [
        Contig(contig_id="c1", sequence="ATCG" * 250, coverage=10.0),
        Contig(contig_id="c2", sequence="GCTA" * 125, coverage=20.0),
        Contig(contig_id="c3", sequence="AAAA" * 75, coverage=5.0),
    ]


@pytest.fixture
def reference_for_profiling() -> tuple[dict[str, str], dict[str, list[tuple[str, str]]]]:
    """Reference sequences and taxonomy for community profiling."""
    ref_seqs = {
        "ref1": "ATCGATCGATCGATCGATCGATCGATCGATCG",
        "ref2": "GGGGCCCCTTTTAAAAGGGCCCCTTTTAAAAG",
    }
    ref_tax = {
        "ref1": [("domain", "Bacteria"), ("phylum", "Firmicutes")],
        "ref2": [("domain", "Bacteria"), ("phylum", "Proteobacteria")],
    }
    return ref_seqs, ref_tax


# ---------------------------------------------------------------------------
# Tests: assemble_contigs
# ---------------------------------------------------------------------------


class TestAssembleContigs:
    """Tests for de Bruijn graph assembly."""

    def test_basic_assembly_from_list(self, overlapping_reads: list[str]) -> None:
        contigs = assemble_contigs(overlapping_reads, k_range=[11], min_contig_length=11)
        assert isinstance(contigs, list)
        # We should get at least one contig
        for c in contigs:
            assert isinstance(c, Contig)
            assert c.length >= 11

    def test_assembly_from_dict(self) -> None:
        reads = {
            "r1": "ATCGATCGATCGATCGATCGATCG",
            "r2": "ATCGATCGATCGATCGATCGATCG",
            "r3": "ATCGATCGATCGATCGATCGATCG",
        }
        contigs = assemble_contigs(reads, k_range=[11], min_contig_length=11, min_kmer_coverage=1)
        assert isinstance(contigs, list)

    def test_empty_reads_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            assemble_contigs([])

    def test_small_k_raises(self) -> None:
        with pytest.raises(ValueError, match="must be >= 11"):
            assemble_contigs(["ATCGATCG"], k_range=[5])

    def test_contig_gc_content(self) -> None:
        reads = ["GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"] * 5
        contigs = assemble_contigs(reads, k_range=[11], min_contig_length=11, min_kmer_coverage=1)
        for c in contigs:
            assert 0.9 <= c.gc_content <= 1.0


# ---------------------------------------------------------------------------
# Tests: calculate_assembly_stats
# ---------------------------------------------------------------------------


class TestCalculateAssemblyStats:
    """Tests for assembly statistics."""

    def test_basic_stats(self, simple_contigs: list[Contig]) -> None:
        stats = calculate_assembly_stats(simple_contigs)
        assert isinstance(stats, AssemblyStats)
        assert stats.total_contigs == 3
        assert stats.total_length == 1000 + 500 + 300
        assert stats.largest_contig == 1000
        assert stats.smallest_contig == 300
        assert stats.n50 > 0
        assert stats.l50 > 0

    def test_n50_single_contig(self) -> None:
        contigs = [Contig("c1", "A" * 1000)]
        stats = calculate_assembly_stats(contigs)
        assert stats.n50 == 1000
        assert stats.l50 == 1

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            calculate_assembly_stats([])

    def test_gc_content_all_gc(self) -> None:
        contigs = [Contig("c1", "GCGCGCGCGCGCGCGCGCGC")]
        stats = calculate_assembly_stats(contigs)
        assert stats.gc_content == pytest.approx(1.0, abs=0.01)


# ---------------------------------------------------------------------------
# Tests: scaffold_contigs
# ---------------------------------------------------------------------------


class TestScaffoldContigs:
    """Tests for scaffolding."""

    def test_basic_scaffolding_no_reads(self, simple_contigs: list[Contig]) -> None:
        scaffolds = scaffold_contigs(simple_contigs)
        assert isinstance(scaffolds, list)
        assert len(scaffolds) == len(simple_contigs)
        for s in scaffolds:
            assert isinstance(s, Scaffold)
            assert s.total_length > 0

    def test_empty_returns_empty(self) -> None:
        assert scaffold_contigs([]) == []


# ---------------------------------------------------------------------------
# Tests: calculate_tetranucleotide_freq
# ---------------------------------------------------------------------------


class TestCalculateTetranucleotideFreq:
    """Tests for TNF calculation."""

    def test_returns_256_values(self) -> None:
        tnf = calculate_tetranucleotide_freq("ATCGATCGATCGATCG")
        assert len(tnf) == 256

    def test_normalized_sums_to_one(self) -> None:
        tnf = calculate_tetranucleotide_freq("ATCGATCGATCGATCG", normalize=True)
        assert abs(sum(tnf) - 1.0) < 0.01

    def test_unnormalized_has_counts(self) -> None:
        tnf = calculate_tetranucleotide_freq("AAAAAAAAAA", normalize=False)
        # AAAA should have high count (both strands counted, so AAAA fwd + TTTT rev)
        assert tnf[0] > 0  # AAAA is index 0

    def test_short_sequence_raises(self) -> None:
        with pytest.raises(ValueError, match="at least 4"):
            calculate_tetranucleotide_freq("ATG")

    def test_homo_polymer_profile(self) -> None:
        tnf_a = calculate_tetranucleotide_freq("AAAAAAAAAA")
        calculate_tetranucleotide_freq("TTTTTTTTTT")
        # Due to reverse complement, poly-A and poly-T should be similar
        assert tnf_a[0] > 0  # AAAA


# ---------------------------------------------------------------------------
# Tests: bin_contigs
# ---------------------------------------------------------------------------


class TestBinContigs:
    """Tests for metagenomic binning."""

    def test_composition_binning(self) -> None:
        contigs = {
            "c1": "ATCG" * 500,
            "c2": "GCTA" * 500,
            "c3": "AAAA" * 500,
        }
        result = bin_contigs(contigs, method="composition", n_bins=2, min_contig_length=100)
        assert isinstance(result, BinningResult)
        assert result.total_contigs == 3
        assert result.binned_contigs > 0
        assert len(result.bins) > 0

    def test_combined_binning_with_coverage(self) -> None:
        contigs = {
            "c1": "ATCG" * 500,
            "c2": "GCTA" * 500,
            "c3": "AAAA" * 500,
        }
        coverage = {"c1": 10.0, "c2": 10.5, "c3": 50.0}
        result = bin_contigs(contigs, coverage, method="combined", n_bins=2, min_contig_length=100)
        assert len(result.bins) == 2

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            bin_contigs({})

    def test_invalid_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown binning"):
            bin_contigs({"c1": "ATCG" * 500}, method="invalid")

    def test_short_contigs_unbinned(self) -> None:
        contigs = {"c1": "ATCG"}  # 4 bp < 1000 default
        result = bin_contigs(contigs, min_contig_length=1000)
        assert len(result.bins) == 0
        assert len(result.unbinned_contigs) == 1


# ---------------------------------------------------------------------------
# Tests: refine_bins
# ---------------------------------------------------------------------------


class TestRefineBins:
    """Tests for bin refinement."""

    def test_basic_refinement(self) -> None:
        contigs = {"c1": "ATCG" * 500, "c2": "GCTA" * 500}
        bins = [
            GenomeBin(
                bin_id="b1",
                contig_ids=["c1", "c2"],
                total_length=4000,
                num_contigs=2,
            )
        ]
        refined = refine_bins(
            bins,
            contigs,
            completeness_threshold=0.0,
            contamination_threshold=1.0,
            marker_genes={"c1": EXPECTED_MARKERS},
            marker_set=EXPECTED_MARKERS,
        )
        assert isinstance(refined, list)
        # At least one bin should survive with these lenient thresholds
        assert len(refined) >= 1

    def test_empty_bins(self) -> None:
        assert refine_bins([], {}) == []


# ---------------------------------------------------------------------------
# Tests: assess_bin_quality
# ---------------------------------------------------------------------------


class TestAssessBinQuality:
    """Tests for bin quality assessment."""

    def test_quality_unavailable_without_marker_genes(self) -> None:
        bins = [GenomeBin(bin_id="b1", contig_ids=["c1"], total_length=3_000_000)]
        assessed = assess_bin_quality(bins)
        assert len(assessed) == 1
        assert assessed[0].completeness is None
        assert assessed[0].contamination is None
        assert assessed[0].quality_score is None

    def test_sequences_alone_do_not_yield_quality(self) -> None:
        contigs = {"c1": "ATCGATCGATCG" * 1000}
        bins = [GenomeBin(bin_id="b1", contig_ids=["c1"], total_length=12000, num_contigs=1)]
        assessed = assess_bin_quality(bins, contigs)
        assert assessed[0].completeness is None
        assert assessed[0].quality_score is None

    def test_empty(self) -> None:
        assert assess_bin_quality([]) == []


# ---------------------------------------------------------------------------
# Tests: build_kmer_index
# ---------------------------------------------------------------------------


class TestBuildKmerIndex:
    """Tests for k-mer index construction."""

    def test_basic_index(self) -> None:
        refs = {"ref1": "ATCGATCGATCGATCGATCGATCGATCGATCG"}
        tax = {"ref1": [("domain", "Bacteria"), ("phylum", "Firmicutes")]}
        index = build_kmer_index(refs, tax, k=11)
        assert isinstance(index, KmerIndex)
        assert index.total_kmers > 0
        assert index.k == 11

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            build_kmer_index({})

    def test_small_k_raises(self) -> None:
        with pytest.raises(ValueError, match="must be >= 11"):
            build_kmer_index({"r1": "ATCGATCG"}, k=5)


# ---------------------------------------------------------------------------
# Tests: profile_community
# ---------------------------------------------------------------------------


class TestProfileCommunity:
    """Tests for community profiling."""

    def test_basic_profiling(
        self,
        reference_for_profiling: tuple[dict[str, str], dict[str, list[tuple[str, str]]]],
    ) -> None:
        ref_seqs, ref_tax = reference_for_profiling
        reads = [ref_seqs["ref1"]]  # Use a read identical to ref1
        profile = profile_community(reads, reference_sequences=ref_seqs, reference_taxonomy=ref_tax, k=11)
        assert isinstance(profile, CommunityProfile)
        assert profile.total_reads == 1

    def test_empty_reads_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            profile_community([])

    def test_no_database_or_refs_raises(self) -> None:
        with pytest.raises(ValueError, match="database or reference"):
            profile_community(["ATCGATCG"])


# ---------------------------------------------------------------------------
# Tests: calculate_relative_abundance
# ---------------------------------------------------------------------------


class TestCalculateRelativeAbundance:
    """Tests for relative abundance normalization."""

    def test_sums_to_one(self) -> None:
        taxa = [
            TaxonProfile("t1", "Firmicutes", "phylum", read_count=70),
            TaxonProfile("t2", "Proteobacteria", "phylum", read_count=30),
        ]
        profile = CommunityProfile(taxa=taxa, classified_reads=100)
        abundances = calculate_relative_abundance(profile)
        assert abs(sum(abundances.values()) - 1.0) < 0.01

    def test_empty_returns_empty(self) -> None:
        profile = CommunityProfile(taxa=[], classified_reads=0)
        assert calculate_relative_abundance(profile) == {}

    def test_min_abundance_filter(self) -> None:
        taxa = [
            TaxonProfile("t1", "Firmicutes", "phylum", read_count=95),
            TaxonProfile("t2", "Rare", "phylum", read_count=5),
        ]
        profile = CommunityProfile(taxa=taxa, classified_reads=100)
        abundances = calculate_relative_abundance(profile, min_abundance=0.10)
        # Rare should be grouped as "Other"
        assert "Firmicutes" in abundances
        assert abs(sum(abundances.values()) - 1.0) < 0.01


class TestTetranucleotideStrandCounting:
    """TNF must count windows on both strands."""

    def test_unnormalized_counts_both_strands(self) -> None:
        # 10 bases -> 7 windows per strand -> total 14 raw counts
        tnf = calculate_tetranucleotide_freq("AAAAAAAAAA", normalize=False)
        assert sum(tnf) == 14.0


class TestRelativeAbundanceRankAggregation:
    """Rank aggregation of relative abundances."""

    def test_aggregates_by_phylum(self) -> None:
        taxa = [
            TaxonProfile(
                "t1",
                "Firmicutes",
                "phylum",
                lineage=[("domain", "Bacteria"), ("phylum", "Firmicutes")],
                read_count=60,
            ),
            TaxonProfile(
                "t2",
                "Bacilli",
                "class",
                lineage=[("domain", "Bacteria"), ("phylum", "Firmicutes"), ("class", "Bacilli")],
                read_count=20,
            ),
            TaxonProfile(
                "t3",
                "Proteobacteria",
                "phylum",
                lineage=[("domain", "Bacteria"), ("phylum", "Proteobacteria")],
                read_count=20,
            ),
        ]
        profile = CommunityProfile(taxa=taxa, classified_reads=100)
        abundances = calculate_relative_abundance(profile, rank="phylum")
        assert abundances["Firmicutes"] == pytest.approx(0.8)
        assert abundances["Proteobacteria"] == pytest.approx(0.2)


class TestScaffoldWithPairedReads:
    """Paired-end links should join contigs into one scaffold."""

    def test_links_join_linked_contigs(self, simple_contigs: list[Contig]) -> None:
        c1, c2, c3 = simple_contigs
        paired_reads = {
            "pair1": (c1.sequence, c2.sequence),
            "pair2": (c1.sequence, c2.sequence),
        }
        scaffolds = scaffold_contigs(simple_contigs, paired_reads=paired_reads, min_links=1)
        assert len(scaffolds) == 2  # {c1, c2} joined; c3 alone
        joined = max(scaffolds, key=lambda s: s.total_length)
        assert len(joined.contigs) == 2
        assert joined.total_length == c1.length + c2.length + joined.gaps[0]


class TestAssemblyStatsDerived:
    """Derived statistics computed from known contigs."""

    def test_mean_coverage_and_median_length(self, simple_contigs: list[Contig]) -> None:
        stats = calculate_assembly_stats(simple_contigs)
        assert stats.mean_coverage == pytest.approx((10.0 + 20.0 + 5.0) / 3.0)
        assert stats.median_length == pytest.approx(500.0)


class TestProfileCommunityUnclassified:
    """Reads with no database support stay unclassified."""

    def test_unrelated_read_unclassified(
        self,
        reference_for_profiling: tuple[dict[str, str], dict[str, list[tuple[str, str]]]],
    ) -> None:
        ref_seqs, ref_tax = reference_for_profiling
        reads = ["TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"]
        profile = profile_community(reads, reference_sequences=ref_seqs, reference_taxonomy=ref_tax, k=11)
        assert profile.classified_reads == 0
        assert profile.unclassified_reads == 1
        assert profile.classification_rate == 0.0
        assert profile.taxa == []


# ---------------------------------------------------------------------------
# Regression tests: marker-gene based bin quality (no fabricated completeness)
# ---------------------------------------------------------------------------


class TestBinQualityFromMarkerGenes:
    """Completeness and contamination come only from supplied marker annotations."""

    def test_three_kb_bin_is_never_scored_complete_without_markers(self) -> None:
        """Any 3 kb bin of ordinary sequence used to score completeness 1.0 from composition alone."""
        contigs = {"c1": _random_dna(3000, random.Random(1))}
        genome_bin = GenomeBin("b1", ["c1"], total_length=3000, num_contigs=1)
        assess_bin_quality([genome_bin], contigs)
        assert genome_bin.completeness is None
        assert genome_bin.contamination is None
        assert genome_bin.quality_score is None
        assert genome_bin.marker_counts == {}

    def test_genome_size_alone_does_not_make_a_bin_complete(self) -> None:
        genome_bin = GenomeBin("b1", ["c1"], total_length=3_500_000)
        assess_bin_quality([genome_bin])
        assert genome_bin.completeness is None and genome_bin.contamination is None

    def test_hand_computed_completeness_contamination_and_score(self) -> None:
        markers = {
            "c1": ["m1", "m2", "m3", "m4", "m5"],
            "c2": ["m6", "m7", "m8", "m9", "m1", "m2", "m2"],  # m1 gains 1 extra copy, m2 gains 2
            "c3": ["m10", "not_a_marker"],
        }
        bin_a = GenomeBin("a", ["c1", "c2"], total_length=3000)
        bin_b = GenomeBin("b", ["c3"], total_length=3000)
        assess_bin_quality([bin_a, bin_b], marker_genes=markers, marker_set=EXPECTED_MARKERS)

        # 9 of 10 markers present; extra copies = (2 - 1) + (3 - 1) = 3
        assert bin_a.completeness == pytest.approx(0.9)
        assert bin_a.contamination == pytest.approx(0.3)
        assert bin_a.quality_score == pytest.approx(0.9 - 5 * 0.3)
        assert bin_a.marker_counts == {"m1": 2, "m2": 3, "m3": 1, "m4": 1, "m5": 1, "m6": 1, "m7": 1, "m8": 1, "m9": 1}
        # markers outside the expected set are ignored
        assert bin_b.completeness == pytest.approx(0.1)
        assert bin_b.contamination == 0.0
        assert bin_b.marker_counts == {"m10": 1}

    def test_complete_bin_with_default_marker_set(self) -> None:
        from metainformant.metagenomics.shotgun.binning import BACTERIAL_SINGLE_COPY_MARKERS

        genome_bin = GenomeBin("b", ["c1"], total_length=3_000_000)
        assess_bin_quality([genome_bin], marker_genes={"c1": list(BACTERIAL_SINGLE_COPY_MARKERS)})
        assert genome_bin.completeness == 1.0
        assert genome_bin.contamination == 0.0

    def test_empty_annotation_means_no_markers_found(self) -> None:
        """An empty mapping is a real result (annotation ran, found nothing), unlike None."""
        genome_bin = GenomeBin("b", ["c1"], total_length=3_000_000)
        assess_bin_quality([genome_bin], marker_genes={}, marker_set=EXPECTED_MARKERS)
        assert genome_bin.completeness == 0.0
        assert genome_bin.contamination == 0.0

    def test_empty_marker_set_raises(self) -> None:
        genome_bin = GenomeBin("b", ["c1"])
        with pytest.raises(ValueError, match="at least one marker"):
            assess_bin_quality([genome_bin], marker_genes={"c1": ["m1"]}, marker_set=[])

    def test_contamination_can_exceed_one(self) -> None:
        """CheckM counts every extra copy, so triplicated markers give 200% contamination."""
        genome_bin = GenomeBin("b", ["c1"])
        assess_bin_quality([genome_bin], marker_genes={"c1": EXPECTED_MARKERS * 3}, marker_set=EXPECTED_MARKERS)
        assert genome_bin.completeness == 1.0
        assert genome_bin.contamination == pytest.approx(2.0)

    def test_reassessment_without_markers_clears_stale_values(self) -> None:
        genome_bin = GenomeBin("b", ["c1"])
        assess_bin_quality([genome_bin], marker_genes={"c1": EXPECTED_MARKERS}, marker_set=EXPECTED_MARKERS)
        assert genome_bin.completeness == 1.0
        assess_bin_quality([genome_bin])
        assert genome_bin.completeness is None and genome_bin.marker_counts == {}


class TestHighQualityBinCount:
    """BinningResult.high_quality_bins used to stay at 0 regardless of the bins."""

    CONTIGS = {"c1": "ATCG" * 500, "c2": "GCTA" * 500, "c3": "AAAA" * 500}

    def test_bin_contigs_fills_high_quality_bins_from_markers(self) -> None:
        markers = {
            "c1": EXPECTED_MARKERS,  # complete, clean: high quality
            "c2": EXPECTED_MARKERS[:6],  # 60% complete
            "c3": EXPECTED_MARKERS + ["m1"],  # complete but 10% contaminated
        }
        result = bin_contigs(
            self.CONTIGS,
            method="composition",
            n_bins=3,
            min_contig_length=100,
            marker_genes=markers,
            marker_set=EXPECTED_MARKERS,
        )
        assert len(result.bins) == 3
        assert result.high_quality_bins == 1
        assert sorted(b.completeness for b in result.bins if b.completeness is not None) == [0.6, 1.0, 1.0]

    def test_high_quality_bins_unavailable_without_markers(self) -> None:
        result = bin_contigs(self.CONTIGS, method="composition", n_bins=3, min_contig_length=100)
        assert result.high_quality_bins is None
        assert all(b.completeness is None and b.contamination is None for b in result.bins)

    def test_no_bins_means_zero_high_quality_bins(self) -> None:
        result = bin_contigs({"c1": "ATCG"}, min_contig_length=1000)
        assert result.high_quality_bins == 0

    def test_count_high_quality_bins_thresholds(self) -> None:
        from metainformant.metagenomics.shotgun.binning import count_high_quality_bins

        bins = [
            GenomeBin("a", completeness=0.95, contamination=0.02),
            GenomeBin("b", completeness=0.90, contamination=0.05),  # boundary values qualify
            GenomeBin("c", completeness=0.89, contamination=0.01),
            GenomeBin("d", completeness=1.0, contamination=0.06),
        ]
        assert count_high_quality_bins(bins) == 2
        assert count_high_quality_bins(bins, min_completeness=0.5, max_contamination=1.0) == 4
        assert count_high_quality_bins([]) == 0

    def test_count_is_unknown_if_any_bin_is_unassessed(self) -> None:
        from metainformant.metagenomics.shotgun.binning import count_high_quality_bins

        assert count_high_quality_bins([GenomeBin("a", completeness=1.0, contamination=0.0), GenomeBin("b")]) is None


class TestRefineBinsWithMarkers:
    """refine_bins needs real marker annotations to filter and split."""

    CONTIGS = {"c1": "ATCG" * 500, "c2": "GCTA" * 500, "c3": "AAAA" * 500}

    def test_requires_marker_genes(self) -> None:
        bins = [GenomeBin("b1", ["c1"], total_length=2000, num_contigs=1)]
        with pytest.raises(ValueError, match="marker_genes"):
            refine_bins(bins, self.CONTIGS)

    def test_drops_incomplete_bins_and_renumbers(self) -> None:
        bins = [
            GenomeBin("a", ["c1"], total_length=2000, num_contigs=1),
            GenomeBin("b", ["c2"], total_length=2000, num_contigs=1),
            GenomeBin("c", ["c3"], total_length=2000, num_contigs=1),
        ]
        markers = {"c1": EXPECTED_MARKERS, "c2": EXPECTED_MARKERS[:3], "c3": EXPECTED_MARKERS[:8]}
        refined = refine_bins(
            bins,
            self.CONTIGS,
            completeness_threshold=0.5,
            contamination_threshold=0.1,
            marker_genes=markers,
            marker_set=EXPECTED_MARKERS,
        )
        assert [b.completeness for b in refined] == [1.0, 0.8]
        assert [b.bin_id for b in refined] == ["refined_bin_0000", "refined_bin_0001"]

    def test_contaminated_bin_is_split_and_subbins_are_assessed(self) -> None:
        contigs = {
            "x1": "ATCG" * 500,
            "x2": "ATCG" * 400,
            "y1": "GGGGCCCCAAAATTTT" * 125,
            "y2": "GGGGCCCCAAAATTTT" * 100,
        }
        # Each composition group carries one complete marker set, so the merged bin has every marker twice
        markers = {
            "x1": EXPECTED_MARKERS[:5],
            "x2": EXPECTED_MARKERS[5:],
            "y1": EXPECTED_MARKERS[:5],
            "y2": EXPECTED_MARKERS[5:],
        }
        merged = GenomeBin("m", list(contigs), total_length=sum(len(v) for v in contigs.values()), num_contigs=4)
        refined = refine_bins(
            [merged],
            contigs,
            completeness_threshold=0.9,
            contamination_threshold=0.1,
            marker_genes=markers,
            marker_set=EXPECTED_MARKERS,
        )
        assert len(refined) == 2
        assert sorted(sorted(b.contig_ids) for b in refined) == [["x1", "x2"], ["y1", "y2"]]
        assert all(b.completeness == 1.0 and b.contamination == 0.0 for b in refined)


# ---------------------------------------------------------------------------
# Regression tests: scaffolding anchors reads on both strands
# ---------------------------------------------------------------------------


def _subprocess_env(**extra: str) -> dict[str, str]:
    """Environment for a child interpreter that imports the same ``metainformant`` as this process.

    ``tests/conftest.py`` prepends the checkout's ``src`` to ``PYTHONPATH``, so inheriting the
    variable could silently select a different copy of the package than the one under test.
    """
    import metainformant

    package_root = os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__)))
    return {**os.environ, "PYTHONPATH": package_root, **extra}


def _reverse_complement(seq: str) -> str:
    """Independent oracle for the reverse complement."""
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


class TestScaffoldAnchorsBothStrands:
    """A mate sequenced from the opposite strand maps to the reverse complement of its contig."""

    @staticmethod
    def _contigs() -> tuple[Contig, Contig]:
        rng = random.Random(17)
        return (
            Contig("c1", _random_dna(400, rng), coverage=10.0),
            Contig("c2", _random_dna(400, rng), coverage=10.0),
        )

    def test_reverse_strand_mate_links_contigs(self) -> None:
        c1, c2 = self._contigs()
        # Standard FR pair: read 1 forward on c1, read 2 is the reverse complement of c2
        pairs = {"p1": (c1.sequence[20:80], _reverse_complement(c2.sequence[100:160]))}
        scaffolds = scaffold_contigs([c1, c2], paired_reads=pairs, min_links=1)
        assert len(scaffolds) == 1
        assert {c.contig_id for c in scaffolds[0].contigs} == {"c1", "c2"}

    def test_first_read_on_reverse_strand_links_contigs(self) -> None:
        c1, c2 = self._contigs()
        pairs = {"p1": (_reverse_complement(c1.sequence[200:260]), c2.sequence[10:70])}
        scaffolds = scaffold_contigs([c1, c2], paired_reads=pairs, min_links=1)
        assert len(scaffolds) == 1

    def test_both_reads_on_reverse_strand_link_contigs(self) -> None:
        c1, c2 = self._contigs()
        pairs = {"p1": (_reverse_complement(c1.sequence[5:65]), _reverse_complement(c2.sequence[300:360]))}
        assert len(scaffold_contigs([c1, c2], paired_reads=pairs, min_links=1)) == 1

    def test_min_links_counts_reverse_strand_pairs(self) -> None:
        c1, c2 = self._contigs()
        pairs = {
            f"p{i}": (c1.sequence[i : i + 60], _reverse_complement(c2.sequence[100 + i : 160 + i])) for i in range(3)
        }
        assert len(scaffold_contigs([c1, c2], paired_reads=pairs, min_links=3)) == 1
        assert len(scaffold_contigs([c1, c2], paired_reads=pairs, min_links=4)) == 2

    def test_unrelated_reads_do_not_link(self) -> None:
        c1, c2 = self._contigs()
        other = _random_dna(60, random.Random(99))
        pairs = {"p1": (c1.sequence[20:80], _reverse_complement(other))}
        assert len(scaffold_contigs([c1, c2], paired_reads=pairs, min_links=1)) == 2


class TestAssemblyIsDeterministic:
    """Contig order and IDs must not depend on the per-process string hash seed."""

    CODE = textwrap.dedent("""
        import random
        from metainformant.metagenomics.shotgun.assembly import assemble_contigs

        rng = random.Random(3)
        genomes = ["".join(rng.choice("ACGT") for _ in range(300)) for _ in range(3)]
        reads = [g[i : i + 60] for g in genomes for i in range(0, 241, 5) for _ in range(3)]
        contigs = assemble_contigs(reads, k_range=[21], min_contig_length=100, min_kmer_coverage=2)
        print([(c.contig_id, c.sequence) for c in contigs])
        """)

    def test_output_is_identical_across_hash_seeds(self) -> None:
        children = [
            subprocess.Popen(
                [sys.executable, "-c", self.CODE],
                env=_subprocess_env(PYTHONHASHSEED=hash_seed),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            for hash_seed in ("0", "1")
        ]
        outputs = set()
        for child in children:
            stdout, stderr = child.communicate(timeout=300)
            assert child.returncode == 0, stderr
            outputs.add(stdout)
        assert len(outputs) == 1


# ---------------------------------------------------------------------------
# Regression tests: reads without a common ancestor are not "classified"
# ---------------------------------------------------------------------------


class TestProfileNoCommonAncestor:
    """A read whose k-mers span taxa with no shared rank used to be counted as classified."""

    @staticmethod
    def _references() -> tuple[dict[str, str], dict[str, list[tuple[str, str]]], str, str, str]:
        rng = random.Random(9)
        shared = _random_dna(60, rng)  # identical segment present in a bacterium and an archaeon
        bac_only = _random_dna(80, rng)
        arc_only = _random_dna(80, rng)
        refs = {"bac": bac_only + shared, "arc": arc_only + shared}
        tax = {
            "bac": [("domain", "Bacteria"), ("phylum", "Firmicutes")],
            "arc": [("domain", "Archaea"), ("phylum", "Euryarchaeota")],
        }
        return refs, tax, bac_only, arc_only, shared

    def test_split_hits_without_common_ancestor_stay_unclassified(self) -> None:
        refs, tax, bac_only, arc_only, _ = self._references()
        read = bac_only[:45] + arc_only[:45]  # half the k-mers support each domain
        profile = profile_community(
            [read], reference_sequences=refs, reference_taxonomy=tax, k=11, confidence_threshold=0.6
        )
        assert profile.classified_reads == 0
        assert profile.unclassified_reads == 1
        assert profile.classification_rate == 0.0
        assert profile.taxa == []

    def test_kmers_shared_across_domains_do_not_classify_a_read(self) -> None:
        refs, tax, _, _, shared = self._references()
        profile = profile_community([shared], reference_sequences=refs, reference_taxonomy=tax, k=11)
        assert profile.classified_reads == 0
        assert profile.unclassified_reads == 1
        assert all(t.name != "unclassified" for t in profile.taxa)

    def test_uninformative_hits_do_not_veto_a_clear_assignment(self) -> None:
        refs, tax, bac_only, _, shared = self._references()
        read = bac_only[:50] + shared[:30]  # bacterial k-mers plus k-mers shared with the archaeon
        profile = profile_community([read], reference_sequences=refs, reference_taxonomy=tax, k=11)
        assert profile.classified_reads == 1
        assert [t.name for t in profile.taxa] == ["Firmicutes"]

    def test_read_from_one_taxon_is_still_classified(self) -> None:
        refs, tax, bac_only, _, _ = self._references()
        profile = profile_community([bac_only[:60]], reference_sequences=refs, reference_taxonomy=tax, k=11)
        assert profile.classified_reads == 1
        assert profile.taxa[0].lineage == [("domain", "Bacteria"), ("phylum", "Firmicutes")]
