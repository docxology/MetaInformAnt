"""Tests for metagenomics amplicon submodule.

Tests OTU clustering, ASV denoising, chimera detection, taxonomy classification.
Uses real implementations -- NO mocking per project policy.
"""

from __future__ import annotations

import math

import pytest

from metainformant.metagenomics.amplicon.asv_denoising import (
    DenoisingResult,
    ErrorModel,
    denoise_sequences,
    estimate_error_rates,
    merge_paired_reads,
)
from metainformant.metagenomics.amplicon.otu_clustering import (
    OTU,
    ClusteringResult,
    calculate_identity,
    cluster_otus,
    filter_chimeras,
)
from metainformant.metagenomics.amplicon.taxonomy import (
    TaxonomyAssignment,
    TaxonomyNode,
    build_taxonomy_tree,
    calculate_confidence,
    classify_taxonomy,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def identical_sequences() -> dict[str, str]:
    """Three identical sequences."""
    seq = "ATCGATCGATCGATCGATCG"
    return {"s1": seq, "s2": seq, "s3": seq}


@pytest.fixture
def dissimilar_sequences() -> dict[str, str]:
    """Sequences with low identity to each other."""
    return {
        "s1": "ATCGATCGATCGATCGATCG",
        "s2": "TTTTTTTTTTTTTTTTTTTT",
        "s3": "CCCCCCCCCCCCCCCCCCCC",
    }


@pytest.fixture
def reference_taxonomy() -> tuple[dict[str, str], dict[str, list[tuple[str, str]]]]:
    """Reference database with taxonomy."""
    ref_db = {
        "ref1": "ATCGATCGATCGATCGATCGATCG",
        "ref2": "TTTTTTTTGGGGGGGGCCCCCCCC",
        "ref3": "AAAACCCCGGGGTTTTAAAACCCC",
    }
    ref_tax = {
        "ref1": [("domain", "Bacteria"), ("phylum", "Firmicutes"), ("class", "Bacilli")],
        "ref2": [("domain", "Bacteria"), ("phylum", "Proteobacteria"), ("class", "Gammaproteobacteria")],
        "ref3": [("domain", "Archaea"), ("phylum", "Euryarchaeota"), ("class", "Methanobacteria")],
    }
    return ref_db, ref_tax


# ---------------------------------------------------------------------------
# Tests: calculate_identity
# ---------------------------------------------------------------------------


class TestCalculateIdentity:
    """Tests for pairwise sequence identity calculation."""

    def test_identical_sequences_return_one(self) -> None:
        """Identical sequences must have identity = 1.0."""
        assert calculate_identity("ATCGATCG", "ATCGATCG") == 1.0

    def test_different_sequences_less_than_one(self) -> None:
        """Different sequences must have identity < 1.0."""
        identity = calculate_identity("ATCGATCG", "TTTTTTTT")
        assert identity < 1.0

    def test_single_mismatch(self) -> None:
        """One mismatch in 8 bases should give identity near 7/8."""
        identity = calculate_identity("ATCGATCG", "ATCGTTCG")
        assert 0.5 < identity < 1.0

    def test_empty_sequence_raises(self) -> None:
        """Empty sequences should raise ValueError."""
        with pytest.raises(ValueError, match="non-empty"):
            calculate_identity("", "ATCG")

    def test_empty_second_sequence_raises(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            calculate_identity("ATCG", "")

    def test_strips_gaps(self) -> None:
        """Gaps should be stripped before comparison."""
        identity = calculate_identity("ATC-GATCG", "ATCGATCG")
        assert identity == 1.0

    def test_case_insensitive(self) -> None:
        """Identity calculation should be case-insensitive."""
        assert calculate_identity("atcgatcg", "ATCGATCG") == 1.0


# ---------------------------------------------------------------------------
# Tests: cluster_otus
# ---------------------------------------------------------------------------


class TestClusterOtus:
    """Tests for OTU clustering."""

    def test_identical_sequences_one_otu(self, identical_sequences: dict[str, str]) -> None:
        """Identical sequences should cluster into one OTU."""
        result = cluster_otus(identical_sequences, threshold=0.97)
        assert isinstance(result, ClusteringResult)
        assert result.num_otus == 1
        assert result.total_sequences == 3

    def test_dissimilar_sequences_multiple_otus(self, dissimilar_sequences: dict[str, str]) -> None:
        """Very dissimilar sequences should form separate OTUs."""
        result = cluster_otus(dissimilar_sequences, threshold=0.97)
        assert result.num_otus > 1
        assert result.total_sequences == 3

    def test_otu_table_populated(self, identical_sequences: dict[str, str]) -> None:
        """OTU table should map centroids to sizes."""
        result = cluster_otus(identical_sequences, threshold=0.97)
        assert len(result.otu_table) == result.num_otus
        assert sum(result.otu_table.values()) == result.total_sequences

    def test_invalid_threshold_raises(self) -> None:
        """Threshold outside (0, 1] should raise ValueError."""
        seqs = {"s1": "ATCGATCG"}
        with pytest.raises(ValueError, match="Threshold"):
            cluster_otus(seqs, threshold=0.0)
        with pytest.raises(ValueError, match="Threshold"):
            cluster_otus(seqs, threshold=1.5)

    def test_empty_sequences_raises(self) -> None:
        """Empty sequences dict should raise ValueError."""
        with pytest.raises(ValueError, match="not be empty"):
            cluster_otus({})

    def test_abundance_sorting(self) -> None:
        """Sorting by abundance should prioritize high-count sequences."""
        seqs = {"rare": "ATCGATCGATCGATCGATCG", "common": "ATCGATCGATCGATCGATCG"}
        abundance = {"rare": 1, "common": 100}
        result = cluster_otus(seqs, threshold=0.97, abundance=abundance, sort_by="abundance")
        assert result.num_otus == 1
        # The centroid should be the common one (sorted first by abundance)
        assert result.otus[0].centroid_id == "common"

    def test_low_threshold_clusters_more(self, dissimilar_sequences: dict[str, str]) -> None:
        """A very low threshold should merge more sequences."""
        high = cluster_otus(dissimilar_sequences, threshold=0.99)
        low = cluster_otus(dissimilar_sequences, threshold=0.50, prefilter=False)
        assert low.num_otus <= high.num_otus


# ---------------------------------------------------------------------------
# Tests: filter_chimeras
# ---------------------------------------------------------------------------


class TestFilterChimeras:
    """Tests for chimera detection."""

    def test_returns_dict_of_bools(self) -> None:
        seqs = {"s1": "AAAAAAAAAA", "s2": "TTTTTTTTTT", "s3": "AAAAATTTTT"}
        result = filter_chimeras(seqs)
        assert isinstance(result, dict)
        assert all(isinstance(v, bool) for v in result.values())
        assert set(result.keys()) == set(seqs.keys())

    def test_most_abundant_not_chimeric(self) -> None:
        """The two most abundant sequences should never be flagged chimeric in de novo mode."""
        seqs = {"s1": "A" * 50, "s2": "T" * 50, "s3": "A" * 25 + "T" * 25}
        abundance = {"s1": 1000, "s2": 500, "s3": 1}
        result = filter_chimeras(seqs, abundance=abundance)
        assert result["s1"] is False
        assert result["s2"] is False

    def test_reference_db_mode(self) -> None:
        """With a reference DB, all sequences are checked against it."""
        seqs = {"q1": "ATCGATCGATCGATCG", "q2": "TTTTTTTTTTTTTTTT"}
        ref_db = {"r1": "ATCGATCGATCGATCG", "r2": "TTTTTTTTTTTTTTTT"}
        result = filter_chimeras(seqs, reference_db=ref_db)
        assert isinstance(result, dict)
        assert len(result) == 2

    def test_empty_sequences_returns_empty(self) -> None:
        assert filter_chimeras({}) == {}


# ---------------------------------------------------------------------------
# Tests: denoise_sequences
# ---------------------------------------------------------------------------


class TestDenoiseSequences:
    """Tests for ASV denoising."""

    def test_basic_denoising(self) -> None:
        reads = {
            "r1": "ATCGATCGATCG",
            "r2": "ATCGATCGATCG",
            "r3": "ATCGATCGATCG",
            "r4": "ATCGATCGATCG",
            "r5": "ATCGATTGATCG",  # one substitution
        }
        result = denoise_sequences(reads)
        assert isinstance(result, DenoisingResult)
        assert result.num_asvs >= 1
        assert result.total_reads >= 4

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            denoise_sequences({})


# ---------------------------------------------------------------------------
# Tests: estimate_error_rates
# ---------------------------------------------------------------------------


class TestEstimateErrorRates:
    """Tests for error rate estimation."""

    def test_basic_operation(self) -> None:
        scores = [[30, 30, 25, 20], [28, 30, 22, 18]]
        model = estimate_error_rates(scores)
        assert isinstance(model, ErrorModel)
        assert len(model.transition_rates) == 4
        assert model.quality_error_rates[30] == pytest.approx(0.001, rel=1e-3)

    def test_with_sequences(self) -> None:
        scores = [[30, 30, 30, 30], [30, 30, 30, 30]]
        seqs = ["ATCG", "ATCG"]
        model = estimate_error_rates(scores, sequences=seqs)
        assert model.positions_modeled == 4

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            estimate_error_rates([])


# ---------------------------------------------------------------------------
# Tests: merge_paired_reads
# ---------------------------------------------------------------------------


class TestMergePairedReads:
    """Tests for paired-end read merging."""

    def test_basic_merge(self) -> None:
        overlap = "ATCGATCGATCG"
        fwd = {"r1": "AAAAAAAAAA" + overlap}
        # reverse complement of overlap region at start of reverse read
        rev = {"r1": overlap + "GGGGGGGGGG"}
        merged = merge_paired_reads(fwd, rev, min_overlap=5)
        # Should be able to merge if the overlap region matches
        assert isinstance(merged, dict)

    def test_no_common_ids_raises(self) -> None:
        with pytest.raises(ValueError, match="No common read IDs"):
            merge_paired_reads({"r1": "ATCG"}, {"r2": "ATCG"})


# ---------------------------------------------------------------------------
# Tests: classify_taxonomy
# ---------------------------------------------------------------------------


class TestClassifyTaxonomy:
    """Tests for taxonomic classification."""

    def test_basic_classification(
        self,
        reference_taxonomy: tuple[dict[str, str], dict[str, list[tuple[str, str]]]],
    ) -> None:
        ref_db, ref_tax = reference_taxonomy
        query = {"q1": "ATCGATCGATCGATCGATCGATCG"}
        results = classify_taxonomy(query, ref_db, ref_tax, method="naive_bayes", bootstrap_n=10)
        assert len(results) == 1
        assert isinstance(results[0], TaxonomyAssignment)
        assert results[0].sequence_id == "q1"

    def test_blast_method(
        self,
        reference_taxonomy: tuple[dict[str, str], dict[str, list[tuple[str, str]]]],
    ) -> None:
        ref_db, ref_tax = reference_taxonomy
        query = {"q1": "ATCGATCGATCGATCGATCGATCG"}
        results = classify_taxonomy(query, ref_db, ref_tax, method="blast")
        assert len(results) == 1
        assert results[0].method == "blast"

    def test_empty_sequences_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            classify_taxonomy({}, {"r1": "ATCG"})

    def test_empty_reference_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            classify_taxonomy({"q1": "ATCG"}, {})

    def test_invalid_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown classification"):
            classify_taxonomy({"q1": "ATCG"}, {"r1": "ATCG"}, method="invalid")


# ---------------------------------------------------------------------------
# Tests: build_taxonomy_tree
# ---------------------------------------------------------------------------


class TestBuildTaxonomyTree:
    """Tests for taxonomy tree construction."""

    def test_basic_tree(self) -> None:
        assignments = [
            TaxonomyAssignment("s1", [("domain", "Bacteria"), ("phylum", "Firmicutes")]),
            TaxonomyAssignment("s2", [("domain", "Bacteria"), ("phylum", "Proteobacteria")]),
        ]
        tree = build_taxonomy_tree(assignments)
        assert isinstance(tree, TaxonomyNode)
        assert tree.rank == "root"
        assert tree.total_descendants() == 2
        assert "Bacteria" in tree.children

    def test_tree_to_dict(self) -> None:
        assignments = [
            TaxonomyAssignment("s1", [("domain", "Bacteria")]),
        ]
        tree = build_taxonomy_tree(assignments)
        d = tree.to_dict()
        assert "name" in d
        assert d["total"] == 1


# ---------------------------------------------------------------------------
# Tests: calculate_confidence
# ---------------------------------------------------------------------------


class TestCalculateConfidence:
    """Tests for confidence calculation."""

    def test_returns_float_stats(self) -> None:
        assignments = [
            TaxonomyAssignment("s1", [], confidence={"domain": 0.95, "phylum": 0.80}),
            TaxonomyAssignment("s2", [], confidence={"domain": 0.90, "phylum": 0.60}),
        ]
        stats = calculate_confidence(assignments, min_confidence=0.8)
        assert "domain" in stats
        assert stats["domain"]["mean"] == pytest.approx(0.925, abs=1e-6)
        assert 0.0 <= stats["domain"]["fraction_confident"] <= 1.0

    def test_empty_assignments_returns_empty(self) -> None:
        assert calculate_confidence([]) == {}
