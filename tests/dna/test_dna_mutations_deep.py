"""Deep tests for metainformant.dna.variation.mutations (real computation, no test doubles)."""

import pytest

from metainformant.dna.variation import mutations


class TestCalculateMutationRate:
    def test_identical_sequences_zero_rate(self) -> None:
        assert mutations.calculate_mutation_rate("ATCGATCG", "ATCGATCG") == 0.0

    def test_known_rate(self) -> None:
        # 2 differences out of 8 sites = 0.25
        assert mutations.calculate_mutation_rate("AAAAAAAA", "AATAACAA") == pytest.approx(0.25)

    def test_case_insensitive(self) -> None:
        # Both sequences uppercased before comparison: "acgt" vs "TCGA" differ at 2/4 sites
        assert mutations.calculate_mutation_rate("acgt", "TCGA") == pytest.approx(0.5)
        assert mutations.calculate_mutation_rate("ACGT", "TCGA") == pytest.approx(0.5)

    def test_empty_sequences(self) -> None:
        assert mutations.calculate_mutation_rate("", "") == 0.0

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="equal length"):
            mutations.calculate_mutation_rate("ATCG", "ATCGAA")


class TestClassifyMutations:
    def test_transition_counted(self) -> None:
        # A -> G is a transition
        result = mutations.classify_mutations("AAAA", "AGAA")
        assert result["total"] == 1
        assert result["transitions"] == 1
        assert result["transversions"] == 0

    def test_transversion_counted(self) -> None:
        # A -> C is a transversion
        result = mutations.classify_mutations("AAAA", "ACAA")
        assert result["total"] == 1
        assert result["transversions"] == 1
        assert result["transitions"] == 0

    def test_no_mutations_all_zero(self) -> None:
        result = mutations.classify_mutations("ATCG", "ATCG")
        assert result["total"] == 0
        assert result["synonymous"] == 0
        assert result["nonsynonymous"] == 0

    def test_synonymous_vs_nonsynonymous_codon_context(self) -> None:
        # Position 2 (0-indexed) of codon TTT -> TTC is synonymous (both F)
        result_syn = mutations.classify_mutations("TTTAAA", "TTCAAA")
        assert result_syn["total"] == 1
        assert result_syn["synonymous"] == 1
        assert result_syn["nonsynonymous"] == 0

        # TTT -> TAT at position 1 changes F -> Y, nonsynonymous
        result_nonsyn = mutations.classify_mutations("TTTAAA", "TATAAA")
        assert result_nonsyn["total"] == 1
        assert result_nonsyn["nonsynonymous"] == 1
        assert result_nonsyn["synonymous"] == 0

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="equal length"):
            mutations.classify_mutations("ATCG", "AT")


class TestFindMutationalHotspots:
    def test_basic_output_structure(self) -> None:
        seq = "ATCGATCGATCG"
        hotspots = mutations.find_mutational_hotspots(seq, window_size=4)
        assert len(hotspots) == len(seq) - 4 + 1
        for pos, score in hotspots:
            assert pos >= 0
            assert 0.0 <= score <= 1.0

    def test_gc_rich_window_scores_higher(self) -> None:
        seq = "GGGGAAAAGGGG"
        hotspots = mutations.find_mutational_hotspots(seq, window_size=4)
        gc_window_score = hotspots[0][1]  # GGGG
        at_window_score = hotspots[2][1]  # AAAA
        assert gc_window_score > at_window_score

    def test_empty_or_invalid_window(self) -> None:
        assert mutations.find_mutational_hotspots("", window_size=10) == []
        assert mutations.find_mutational_hotspots("ATCG", window_size=0) == []


class TestCalculateSubstitutionMatrix:
    def test_counts_substitutions(self) -> None:
        result = mutations.calculate_substitution_matrix("ATCG", "AGCC")
        assert result == {("T", "G"): 1, ("G", "C"): 1}

    def test_identical_sequences_empty_matrix(self) -> None:
        assert mutations.calculate_substitution_matrix("ATCG", "ATCG") == {}

    def test_length_mismatch_returns_empty(self) -> None:
        assert mutations.calculate_substitution_matrix("ATCG", "ATCGAA") == {}


class TestDetectSelectionSignatures:
    def test_identical_to_refs_full_conservation(self) -> None:
        result = mutations.detect_selection_signatures("ATCGATCG", ["ATCGATCG", "ATCGATCG"])
        assert result["diversity"] == pytest.approx(0.0)
        assert result["conservation"] == pytest.approx(1.0)

    def test_empty_refs_default(self) -> None:
        result = mutations.detect_selection_signatures("ATCG", [])
        assert result == {"diversity": 0.0, "conservation": 1.0}

    def test_diversity_between_zero_and_one(self) -> None:
        result = mutations.detect_selection_signatures("ATCGATCG", ["ATCGATCA", "TCGAATCG"])
        assert 0.0 < result["diversity"] < 1.0
        assert result["conservation"] == pytest.approx(1.0 - result["diversity"])

    def test_length_mismatched_refs_ignored(self) -> None:
        # refs of wrong length contribute no diversity samples
        result = mutations.detect_selection_signatures("ATCG", ["ATCGAAA", "TCG"])
        assert result["diversity"] == 0.0


class TestGenerateMutantLibrary:
    def test_basic_library(self) -> None:
        mutants = mutations.generate_mutant_library("ATCG", [1, 3], ["G", "T"])
        assert mutants == ["AGCG", "ATCT"]

    def test_lowercase_mutation_uppercased(self) -> None:
        mutants = mutations.generate_mutant_library("ATCG", [0], ["c"])
        assert mutants == ["CTCG"]

    def test_unequal_lengths_raise(self) -> None:
        with pytest.raises(ValueError, match="equal length"):
            mutations.generate_mutant_library("ATCG", [1], [])

    def test_out_of_range_position_skipped(self) -> None:
        mutants = mutations.generate_mutant_library("ATCG", [1, 99], ["G", "A"])
        assert mutants == ["AGCG"]


class TestAnalyzeMutationSpectrum:
    def test_spectrum_structure(self) -> None:
        spectrum = mutations.analyze_mutation_spectrum("AGCGATCG", "ATCGATCG")
        assert "mutation_types" in spectrum
        assert "mutation_density" in spectrum
        assert "ti_tv_ratio" in spectrum
        assert "sequence_length" in spectrum
        assert "substitution_matrix" in spectrum
        assert spectrum["sequence_length"] == 8
        assert spectrum["mutation_types"]["total"] == 1

    def test_density_and_ratio(self) -> None:
        # A->C (x2, transversion), T->G (x2, transversion), C->A... 6 total diffs
        # (classify_mutations counts position 1 T->G twice via reverse-complement-like sym; measured: 6)
        spectrum = mutations.analyze_mutation_spectrum("ACGCACGC", "ATCGATCG")
        assert spectrum["mutation_density"] == pytest.approx(0.75)
        assert spectrum["mutation_types"]["transitions"] == 2
        assert spectrum["ti_tv_ratio"] == pytest.approx(2 / 4)

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="equal length"):
            mutations.analyze_mutation_spectrum("ATCG", "ATCGAA")


class TestSimulateSequenceEvolution:
    def test_lineage_length(self) -> None:
        lineage = mutations.simulate_sequence_evolution("ATCG", generations=3, mutation_rate=0.9)
        assert len(lineage) == 4  # original + 3 generations
        assert lineage[0] == "ATCG"

    def test_zero_or_negative_generations_returns_original(self) -> None:
        assert mutations.simulate_sequence_evolution("ATCG", generations=0) == ["ATCG"]
        assert mutations.simulate_sequence_evolution("ATCG", generations=-1) == ["ATCG"]

    def test_empty_sequence(self) -> None:
        assert mutations.simulate_sequence_evolution("", generations=5) == [""]

    def test_high_rate_eventually_changes(self) -> None:
        lineage = mutations.simulate_sequence_evolution("AAAAAAAAAA", generations=10, mutation_rate=1.0)
        # With rate 1.0 and num_mutations=1 per generation, later generations differ from start
        assert any(seq != lineage[0] for seq in lineage[1:])
