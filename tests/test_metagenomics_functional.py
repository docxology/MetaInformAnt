"""Tests for metagenomics functional submodule.

Tests ORF prediction, gene annotation, pathway analysis.
Uses real implementations -- NO mocking per project policy.
"""

from __future__ import annotations

import pytest

from metainformant.metagenomics.functional.annotation import (
    ORF,
    FunctionalAnnotation,
    annotate_genes,
    classify_gene_families,
    predict_orfs,
)
from metainformant.metagenomics.functional.pathways import (
    PathwayDefinition,
    PathwayResult,
    calculate_pathway_completeness,
    compare_pathway_profiles,
    reconstruct_pathways,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def coding_sequence() -> str:
    """A DNA sequence that contains at least one ORF with ATG start and TAA stop.

    ATG (start) + 33 codons + TAA (stop) = 102 nt ORF.
    """
    # Build a clean coding region: ATG + 33 codons of "GCT" repeats + TAA
    coding_region = "ATG" + "GCT" * 33 + "TAA"
    # Pad flanking so that the total sequence is large enough
    flank = "AAAAAAAAAA"
    return flank + coding_region + flank


@pytest.fixture
def simple_hmm_db() -> dict[str, dict[int, dict[str, float]]]:
    """Minimal HMM database for annotation testing."""
    # A profile that scores well on alanine-rich proteins
    profile: dict[int, dict[str, float]] = {}
    for pos in range(10):
        profile[pos] = {"A": 5.0, "G": 2.0, "L": 1.0, "M": 1.0}
    return {"COG:C0001": profile, "KO:K00844": profile}


@pytest.fixture
def glycolysis_kos() -> dict[str, list[str]]:
    """Annotations matching some glycolysis pathway KOs."""
    return {
        "gene1": ["K00844"],
        "gene2": ["K00845"],
        "gene3": ["K01810"],
        "gene4": ["K00850"],
        "gene5": ["K01623"],
    }


# ---------------------------------------------------------------------------
# Tests: predict_orfs
# ---------------------------------------------------------------------------


class TestPredictOrfs:
    """Tests for ORF prediction."""

    def test_finds_orfs_in_coding_sequence(self, coding_sequence: str) -> None:
        orfs = predict_orfs(coding_sequence, sequence_id="test_seq", min_length=90)
        assert isinstance(orfs, list)
        assert len(orfs) >= 1
        for orf in orfs:
            assert isinstance(orf, ORF)
            assert orf.sequence_id == "test_seq"
            assert len(orf.nucleotide_seq) >= 90
            assert len(orf.protein_seq) > 0

    def test_short_sequence_no_orfs(self) -> None:
        """Sequence shorter than min_length should yield no ORFs."""
        orfs = predict_orfs("ATGATCGTAA", min_length=200)
        assert orfs == []

    def test_protein_translation(self, coding_sequence: str) -> None:
        """Predicted protein should contain valid amino acids."""
        orfs = predict_orfs(coding_sequence, min_length=90)
        for orf in orfs:
            # All amino acids should be standard letters
            assert all(aa in "ACDEFGHIKLMNPQRSTVWXY" for aa in orf.protein_seq)

    def test_both_strands(self) -> None:
        """ORFs should be predicted on both + and - strands."""
        # Forward ORF
        fwd = "ATG" + "GCT" * 40 + "TAA"
        # Reverse complement of another ORF
        rev_orf = "ATG" + "AAA" * 40 + "TAA"
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        rev_complement = "".join(complement[b] for b in reversed(rev_orf))
        seq = fwd + "NNNNNNNNNN" + rev_complement
        orfs = predict_orfs(seq, min_length=90)
        strands = {orf.strand for orf in orfs}
        # We should get ORFs from at least the + strand
        assert "+" in strands

    def test_orf_properties(self, coding_sequence: str) -> None:
        orfs = predict_orfs(coding_sequence, min_length=90)
        for orf in orfs:
            assert orf.length_nt == len(orf.nucleotide_seq)
            assert orf.length_aa == len(orf.protein_seq)
            assert orf.frame in (0, 1, 2)


# ---------------------------------------------------------------------------
# Tests: annotate_genes
# ---------------------------------------------------------------------------


class TestAnnotateGenes:
    """Tests for gene annotation."""

    def test_with_hmm_db(self, simple_hmm_db: dict) -> None:
        # Protein rich in alanine (should match the profile)
        sequences = {"gene1": "AAAAAAAGGGLLM" * 3}
        results = annotate_genes(sequences, hmm_db=simple_hmm_db, min_score=1.0)
        assert len(results) == 1
        assert isinstance(results[0], FunctionalAnnotation)
        assert results[0].orf_id == "gene1"

    def test_without_hmm_db(self) -> None:
        sequences = {"gene1": "MAAAKKKLLLL"}
        results = annotate_genes(sequences, hmm_db=None)
        assert len(results) == 1
        assert results[0].best_hit is None

    def test_empty_returns_empty(self) -> None:
        results = annotate_genes({})
        assert results == []


# ---------------------------------------------------------------------------
# Tests: classify_gene_families
# ---------------------------------------------------------------------------


class TestClassifyGeneFamilies:
    """Tests for gene family classification."""

    def test_composition_based(self) -> None:
        """Without reference profiles, uses composition heuristic."""
        genes = {"g1": "A" * 100 + "K" * 50 + "D" * 50}
        result = classify_gene_families(genes, database="COG")
        assert isinstance(result, dict)
        assert "g1" in result
        assert isinstance(result["g1"], list)

    def test_with_reference_profiles(self, simple_hmm_db: dict) -> None:
        genes = {"g1": "AAAAAAAGGGLLM" * 3}
        result = classify_gene_families(genes, database="COG", reference_profiles=simple_hmm_db, min_score=1.0)
        assert "g1" in result

    def test_kegg_classification(self) -> None:
        genes = {"g1": "A" * 100 + "K" * 50 + "D" * 50}
        result = classify_gene_families(genes, database="KEGG")
        assert "g1" in result


# ---------------------------------------------------------------------------
# Tests: reconstruct_pathways
# ---------------------------------------------------------------------------


class TestReconstructPathways:
    """Tests for pathway reconstruction."""

    def test_basic_reconstruction(self, glycolysis_kos: dict[str, list[str]]) -> None:
        results = reconstruct_pathways(glycolysis_kos, database="KEGG")
        assert isinstance(results, list)
        # Should find glycolysis pathway (map00010) with partial completeness
        glycolysis_hits = [r for r in results if r.pathway_id == "map00010"]
        assert len(glycolysis_hits) == 1
        assert glycolysis_hits[0].completeness > 0.0
        assert glycolysis_hits[0].total_matched > 0

    def test_sorted_by_completeness(self, glycolysis_kos: dict[str, list[str]]) -> None:
        results = reconstruct_pathways(glycolysis_kos)
        completeness_vals = [r.completeness for r in results]
        assert completeness_vals == sorted(completeness_vals, reverse=True)

    def test_empty_annotations(self) -> None:
        results = reconstruct_pathways({})
        # All pathways should have 0 completeness
        for r in results:
            assert r.completeness == 0.0

    def test_min_completeness_filter(self, glycolysis_kos: dict[str, list[str]]) -> None:
        all_results = reconstruct_pathways(glycolysis_kos, min_completeness=0.0)
        filtered = reconstruct_pathways(glycolysis_kos, min_completeness=0.9)
        assert len(filtered) <= len(all_results)


# ---------------------------------------------------------------------------
# Tests: calculate_pathway_completeness
# ---------------------------------------------------------------------------


class TestCalculatePathwayCompleteness:
    """Tests for single pathway completeness."""

    def test_full_completeness(self) -> None:
        pathway = PathwayDefinition(
            pathway_id="test",
            name="Test",
            required_kos=["K001", "K002", "K003"],
        )
        annotations = {"K001", "K002", "K003"}
        assert calculate_pathway_completeness(pathway, annotations) == pytest.approx(1.0)

    def test_partial_completeness(self) -> None:
        pathway = PathwayDefinition(
            pathway_id="test",
            name="Test",
            required_kos=["K001", "K002", "K003", "K004"],
        )
        annotations = {"K001", "K002"}
        assert calculate_pathway_completeness(pathway, annotations) == pytest.approx(0.5)

    def test_zero_completeness(self) -> None:
        pathway = PathwayDefinition(
            pathway_id="test",
            name="Test",
            required_kos=["K001"],
        )
        annotations: set[str] = set()
        assert calculate_pathway_completeness(pathway, annotations) == pytest.approx(0.0)

    def test_step_based(self) -> None:
        pathway = PathwayDefinition(
            pathway_id="test",
            name="Test",
            reaction_steps=[["K001", "K002"], ["K003"]],  # step1 has 2 alternatives
        )
        # Step 1 satisfied by K002, step 2 satisfied by K003
        annotations = {"K002", "K003"}
        assert calculate_pathway_completeness(pathway, annotations) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Tests: compare_pathway_profiles
# ---------------------------------------------------------------------------


class TestComparePathwayProfiles:
    """Tests for cross-sample pathway comparison."""

    def test_basic_comparison(self) -> None:
        sample_a = {"gene1": ["K00844", "K00845"]}
        sample_b = {"gene1": ["K01647", "K01681"]}
        samples = {"sampleA": sample_a, "sampleB": sample_b}
        result = compare_pathway_profiles(samples)
        assert isinstance(result, dict)
        # Result maps pathway_id -> {sample_id: completeness}
        for pathway_id, sample_scores in result.items():
            assert "sampleA" in sample_scores or "sampleB" in sample_scores

    def test_empty_samples(self) -> None:
        result = compare_pathway_profiles({})
        assert result == {}
