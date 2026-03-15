"""Tests for metagenomics shotgun submodule.

Tests assembly, binning, community profiling.
Uses real implementations -- NO mocking per project policy.
"""

from __future__ import annotations

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
        tnf_t = calculate_tetranucleotide_freq("TTTTTTTTTT")
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
        refined = refine_bins(bins, contigs, completeness_threshold=0.0, contamination_threshold=1.0)
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

    def test_returns_completeness_contamination(self) -> None:
        bins = [GenomeBin(bin_id="b1", contig_ids=["c1"], total_length=3_000_000)]
        assessed = assess_bin_quality(bins)
        assert len(assessed) == 1
        assert 0.0 <= assessed[0].completeness <= 1.0
        assert 0.0 <= assessed[0].contamination <= 1.0

    def test_with_contig_sequences(self) -> None:
        contigs = {"c1": "ATCGATCGATCG" * 1000}
        bins = [GenomeBin(bin_id="b1", contig_ids=["c1"], total_length=12000, num_contigs=1)]
        assessed = assess_bin_quality(bins, contigs)
        assert assessed[0].quality_score == assessed[0].completeness - 5.0 * assessed[0].contamination

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
