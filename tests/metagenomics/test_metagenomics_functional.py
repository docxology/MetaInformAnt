"""Tests for metagenomics functional submodule.

Tests ORF prediction, gene annotation, pathway analysis.
Uses real implementations -- real-implementation policy.
"""

from __future__ import annotations

import math
import random

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


class TestPredictOrfsPartial:
    """Partial ORFs at sequence edges must be flagged."""

    def test_orf_without_stop_is_partial(self) -> None:
        # Start codon + 40 codons, sequence ends before any stop codon
        seq = "ATG" + "GCT" * 40
        orfs = predict_orfs(seq, min_length=90)
        assert len(orfs) >= 1
        assert any(orf.partial for orf in orfs)
        for orf in orfs:
            if orf.partial:
                assert len(orf.nucleotide_seq) >= 90


class TestAnnotateGenesECExtraction:
    """EC numbers embedded in family IDs must be extracted."""

    def test_ec_number_extracted(self, simple_hmm_db: dict) -> None:
        hmm_db = {"FAM;EC:1.1.1.1": simple_hmm_db["COG:C0001"]}
        sequences = {"gene1": "AAAAAAAGGGLLM" * 3}
        results = annotate_genes(sequences, hmm_db=hmm_db, min_score=1.0)
        assert results[0].gene_families == ["FAM;EC:1.1.1.1"]
        assert results[0].ec_numbers == ["1.1.1.1"]


# ---------------------------------------------------------------------------
# Regression tests: ORF prediction (duplicates, translation tables)
# ---------------------------------------------------------------------------


class TestPredictOrfsReportsEachOrfOnce:
    """A start codon at position 0 was also registered as the implicit partial start."""

    @staticmethod
    def _keys(orfs: list[ORF]) -> list[tuple[str, int, int, int]]:
        return [(o.strand, o.start, o.end, o.frame) for o in orfs]

    def test_start_codon_at_position_zero_is_not_duplicated(self) -> None:
        seq = "ATG" + "GCT" * 40 + "TAA"
        orfs = predict_orfs(seq, min_length=90)
        keys = self._keys(orfs)
        assert len(keys) == len(set(keys))
        forward_at_zero = [o for o in orfs if o.strand == "+" and o.start == 0]
        assert len(forward_at_zero) == 1
        assert forward_at_zero[0].end == len(seq)

    def test_partial_rule_matches_explicit_start_when_sequence_begins_with_atg(self) -> None:
        seq = "ATG" + "GCT" * 40 + "TAA"
        with_partial = [o for o in predict_orfs(seq, min_length=90, allow_partial=True) if o.strand == "+"]
        without_partial = [o for o in predict_orfs(seq, min_length=90, allow_partial=False) if o.strand == "+"]
        assert self._keys(with_partial) == self._keys(without_partial)

    def test_unfinished_orf_at_the_end_is_not_duplicated(self) -> None:
        seq = "ATG" + "GCT" * 40  # starts with ATG at 0 and never reaches a stop codon
        keys = self._keys(predict_orfs(seq, min_length=90))
        assert len(keys) == len(set(keys))

    def test_no_duplicates_in_random_sequences(self) -> None:
        rng = random.Random(6)
        for _ in range(20):
            seq = "ATG" + "".join(rng.choice("ACGT") for _ in range(400))
            keys = self._keys(predict_orfs(seq, min_length=30))
            assert len(keys) == len(set(keys))


class TestTranslationTable:
    """`translation_table` used to be ignored (always the standard stop codons)."""

    # ATG, 20 Ala, TGA, 20 Ala, TAA: TGA is a stop in table 11 but tryptophan in table 4
    SEQ = "ATG" + "GCT" * 20 + "TGA" + "GCT" * 20 + "TAA"

    @staticmethod
    def _forward_from_zero(orfs: list[ORF]) -> ORF:
        return next(o for o in orfs if o.strand == "+" and o.start == 0)

    def test_table_11_stops_at_tga(self) -> None:
        orf = self._forward_from_zero(predict_orfs(self.SEQ, min_length=60, translation_table=11))
        assert orf.end == 66
        assert orf.protein_seq == "M" + "A" * 20

    def test_table_4_reads_tga_as_tryptophan(self) -> None:
        orf = self._forward_from_zero(predict_orfs(self.SEQ, min_length=60, translation_table=4))
        assert orf.end == len(self.SEQ)
        assert orf.protein_seq == "M" + "A" * 20 + "W" + "A" * 20

    def test_vertebrate_mitochondrial_agg_is_a_stop(self) -> None:
        seq = "ATG" + "GCT" * 20 + "AGG" + "GCT" * 20 + "TAA"
        in_table_11 = self._forward_from_zero(predict_orfs(seq, min_length=60, translation_table=11))
        in_table_2 = self._forward_from_zero(predict_orfs(seq, min_length=60, translation_table=2))
        assert in_table_11.protein_seq == "M" + "A" * 20 + "R" + "A" * 20  # AGG is arginine
        assert in_table_2.protein_seq == "M" + "A" * 20  # AGG is a stop

    def test_tables_1_and_11_share_amino_acid_assignments(self) -> None:
        proteins = {
            table: self._forward_from_zero(predict_orfs(self.SEQ, min_length=60, translation_table=table)).protein_seq
            for table in (1, 11)
        }
        assert proteins[1] == proteins[11]

    @pytest.mark.parametrize("table", [0, 99, -1, "11"])
    def test_unknown_table_raises(self, table: object) -> None:
        with pytest.raises(ValueError, match="translation_table"):
            predict_orfs(self.SEQ, translation_table=table)  # type: ignore[arg-type]

    def test_unknown_table_raises_even_for_short_sequences(self) -> None:
        with pytest.raises(ValueError, match="translation_table"):
            predict_orfs("ATG", min_length=200, translation_table=99)


# ---------------------------------------------------------------------------
# Regression tests: annotation reports only what it computed
# ---------------------------------------------------------------------------


class TestAnnotateGenesDoesNotFabricateEValues:
    """`e_value_threshold` was ignored and every hit claimed E = 0.0."""

    SEQUENCES = {"gene1": "AAAAAAAGGGLLM" * 3}

    def test_numeric_e_value_threshold_is_rejected(self, simple_hmm_db: dict) -> None:
        with pytest.raises(ValueError, match="e_value_threshold"):
            annotate_genes(self.SEQUENCES, hmm_db=simple_hmm_db, e_value_threshold=1e-5, min_score=1.0)

    def test_hits_do_not_claim_a_perfect_e_value(self, simple_hmm_db: dict) -> None:
        result = annotate_genes(self.SEQUENCES, hmm_db=simple_hmm_db, min_score=1.0)[0]
        assert result.best_hit is not None
        assert math.isnan(result.best_hit.e_value)
        assert math.isnan(result.best_hit.domain_e_value)
        assert math.isnan(result.best_hit.bias)
        assert result.best_hit.score > 1.0

    def test_min_score_still_filters(self, simple_hmm_db: dict) -> None:
        result = annotate_genes(self.SEQUENCES, hmm_db=simple_hmm_db, min_score=50.0)[0]
        assert result.best_hit is None
        assert result.cog_categories == [] and result.kegg_orthologs == []


class TestClassifyGeneFamiliesDoesNotInventLabels:
    """Amino-acid composition cannot identify a COG/KEGG/Pfam family."""

    CHARGED = "A" * 100 + "K" * 50 + "D" * 50
    CYS_GLY = "C" * 20 + "G" * 40 + "A" * 140

    @pytest.mark.parametrize("database", ["COG", "KEGG", "Pfam"])
    def test_no_profiles_means_unannotated(self, database: str) -> None:
        genes = {"charged": self.CHARGED, "cys_gly": self.CYS_GLY, "plain": "MKV" * 30}
        result = classify_gene_families(genes, database=database)
        assert result == {"charged": [], "cys_gly": [], "plain": []}

    def test_profiles_still_assign_families(self, simple_hmm_db: dict) -> None:
        result = classify_gene_families(
            {"g1": "AAAAAAAGGGLLM" * 3}, database="COG", reference_profiles=simple_hmm_db, min_score=1.0
        )
        assert result == {"g1": ["COG:C0001"]}


# ---------------------------------------------------------------------------
# Regression tests: pathway completeness never exceeds 1
# ---------------------------------------------------------------------------


class TestReconstructPathwaysStepCompleteness:
    """Completeness divided matched KOs by the number of steps (6 KOs over 2 steps gave 3.0)."""

    KOS = [f"K0000{i}" for i in range(1, 7)]

    def _definition(self) -> PathwayDefinition:
        return PathwayDefinition(
            pathway_id="pw",
            name="Six KOs over two steps",
            required_kos=self.KOS,
            reaction_steps=[self.KOS[:3], self.KOS[3:]],
        )

    def _run(self, annotated: list[str]):
        results = reconstruct_pathways({"g": annotated}, pathway_definitions={"pw": self._definition()})
        assert len(results) == 1
        return results[0]

    def test_all_kos_present_is_exactly_complete(self) -> None:
        result = self._run(self.KOS)
        assert result.completeness == 1.0
        assert (result.total_matched, result.total_required) == (2, 2)

    def test_one_step_present_is_half_complete(self) -> None:
        result = self._run(self.KOS[:3])
        assert result.completeness == 0.5
        assert (result.total_matched, result.total_required) == (1, 2)
        assert result.missing_kos == self.KOS[3:]

    def test_alternatives_of_one_step_count_once(self) -> None:
        assert self._run(["K00001", "K00002", "K00003"]).completeness == 0.5

    def test_nothing_present(self) -> None:
        result = self._run(["K99999"])
        assert result.completeness == 0.0
        assert result.matched_kos == []

    def test_steps_only_definition_is_scored(self) -> None:
        definition = PathwayDefinition("so", "Steps only", reaction_steps=[["K00001", "K00002"], ["K00003"]])
        result = reconstruct_pathways({"g": ["K00002", "K00003"]}, pathway_definitions={"so": definition})[0]
        assert result.completeness == 1.0
        assert result.matched_kos == ["K00002", "K00003"]

    def test_ec_alternatives_in_steps(self) -> None:
        definition = PathwayDefinition("ec", "EC steps", reaction_steps=[["1.1.1.1", "K00001"], ["2.7.1.1"]])
        result = reconstruct_pathways({"g": ["1.1.1.1"]}, pathway_definitions={"ec": definition})[0]
        assert result.completeness == 0.5
        assert result.matched_ecs == ["1.1.1.1"]
        assert result.missing_ecs == ["2.7.1.1"]

    def test_agrees_with_calculate_pathway_completeness(self) -> None:
        rng = random.Random(12)
        definition = self._definition()
        for _ in range(60):
            annotated = rng.sample(self.KOS, rng.randint(0, 6)) + ["K99999"]
            result = self._run(annotated)
            assert 0.0 <= result.completeness <= 1.0
            assert result.completeness == pytest.approx(calculate_pathway_completeness(definition, set(annotated)))

    def test_cross_sample_profiles_stay_within_unit_interval(self) -> None:
        samples = {"a": {"g": self.KOS}, "b": {"g": self.KOS[:2]}, "c": {"g": []}}
        profile = compare_pathway_profiles(samples, pathway_definitions={"pw": self._definition()})
        assert profile["pw"] == {"a": 1.0, "b": 0.5, "c": 0.0}
