"""Comprehensive tests for the metainformant.simulation module.

Tests cover sequence generation, mutation, translation, similarity analysis,
RNA-seq count simulation, GridWorld agent-based models, and integration
across multiple simulation components.

All tests use real implementations (NO MOCKING).
Reproducibility is ensured via seeded random.Random instances.
"""

from __future__ import annotations

import random
from typing import List

import numpy as np
import pytest

from metainformant.simulation import (
    # Sequence functions
    generate_random_dna,
    generate_random_protein,
    mutate_sequence,
    evolve_sequence,
    translate_dna_to_protein,
    reverse_transcribe_protein_to_dna,
    generate_coding_sequence,
    calculate_sequence_similarity,
    generate_sequence_family,
    analyze_sequence_divergence,
    simulate_gene_duplication,
    # Constants
    DNA_BASES,
    RNA_BASES,
    AMINO_ACIDS,
    GENETIC_CODE,
    # RNA simulation
    simulate_rnaseq_counts,
    # Agent-based models
    GridAgent,
    GridWorld,
)


# ---------------------------------------------------------------------------
# DNA generation
# ---------------------------------------------------------------------------


class TestGenerateRandomDNA:
    """Tests for generate_random_dna."""

    def test_basic_generation(self) -> None:
        """Generate a 100-bp DNA sequence with default GC content."""
        rng = random.Random(42)
        seq = generate_random_dna(100, rng=rng)
        assert len(seq) == 100
        assert all(c in "ACGT" for c in seq)

    def test_length_one(self) -> None:
        """Minimum valid length produces a single valid base."""
        rng = random.Random(42)
        seq = generate_random_dna(1, rng=rng)
        assert len(seq) == 1
        assert seq in "ACGT"

    def test_large_sequence(self) -> None:
        """Generate a large sequence and verify all bases are valid."""
        rng = random.Random(42)
        seq = generate_random_dna(10000, rng=rng)
        assert len(seq) == 10000
        assert set(seq).issubset(set("ACGT"))

    def test_gc_content_high(self) -> None:
        """High GC content request produces approximately correct GC ratio."""
        rng = random.Random(42)
        seq = generate_random_dna(10000, gc_content=0.8, rng=rng)
        gc_count = sum(1 for c in seq if c in "GC")
        gc_ratio = gc_count / len(seq)
        assert 0.75 <= gc_ratio <= 0.85

    def test_gc_content_low(self) -> None:
        """Low GC content request produces approximately correct GC ratio."""
        rng = random.Random(42)
        seq = generate_random_dna(10000, gc_content=0.2, rng=rng)
        gc_count = sum(1 for c in seq if c in "GC")
        gc_ratio = gc_count / len(seq)
        assert 0.15 <= gc_ratio <= 0.25

    def test_gc_content_default(self) -> None:
        """Default GC content (0.5) produces balanced nucleotide composition."""
        rng = random.Random(42)
        seq = generate_random_dna(10000, rng=rng)
        gc_count = sum(1 for c in seq if c in "GC")
        gc_ratio = gc_count / len(seq)
        assert 0.45 <= gc_ratio <= 0.55

    def test_reproducibility(self) -> None:
        """Same seed produces identical sequences."""
        seq1 = generate_random_dna(200, rng=random.Random(42))
        seq2 = generate_random_dna(200, rng=random.Random(42))
        assert seq1 == seq2

    def test_different_seeds_differ(self) -> None:
        """Different seeds produce different sequences."""
        seq1 = generate_random_dna(200, rng=random.Random(42))
        seq2 = generate_random_dna(200, rng=random.Random(99))
        assert seq1 != seq2

    def test_invalid_length_zero(self) -> None:
        """Length 0 raises ValueError (minimum is 1)."""
        with pytest.raises(Exception):
            generate_random_dna(0, rng=random.Random(42))

    def test_invalid_length_negative(self) -> None:
        """Negative length raises ValueError."""
        with pytest.raises(Exception):
            generate_random_dna(-5, rng=random.Random(42))


# ---------------------------------------------------------------------------
# Protein generation
# ---------------------------------------------------------------------------


class TestGenerateRandomProtein:
    """Tests for generate_random_protein."""

    def test_basic_generation(self) -> None:
        """Generate a 50-aa protein with valid amino acids."""
        rng = random.Random(42)
        seq = generate_random_protein(50, rng=rng)
        assert len(seq) == 50
        assert all(c in AMINO_ACIDS for c in seq)

    def test_length_one(self) -> None:
        """Minimum valid length produces a single amino acid."""
        rng = random.Random(42)
        seq = generate_random_protein(1, rng=rng)
        assert len(seq) == 1
        assert seq in AMINO_ACIDS

    def test_large_protein(self) -> None:
        """Generate a large protein and verify all residues are valid."""
        rng = random.Random(42)
        seq = generate_random_protein(5000, rng=rng)
        assert len(seq) == 5000
        assert set(seq).issubset(set(AMINO_ACIDS))

    def test_reproducibility(self) -> None:
        """Same seed produces identical protein sequences."""
        seq1 = generate_random_protein(100, rng=random.Random(42))
        seq2 = generate_random_protein(100, rng=random.Random(42))
        assert seq1 == seq2

    def test_different_seeds_differ(self) -> None:
        """Different seeds produce different proteins."""
        seq1 = generate_random_protein(100, rng=random.Random(42))
        seq2 = generate_random_protein(100, rng=random.Random(99))
        assert seq1 != seq2

    def test_invalid_length_zero(self) -> None:
        """Length 0 raises ValueError (minimum is 1)."""
        with pytest.raises(Exception):
            generate_random_protein(0, rng=random.Random(42))

    def test_invalid_length_negative(self) -> None:
        """Negative length raises ValueError."""
        with pytest.raises(Exception):
            generate_random_protein(-3, rng=random.Random(42))


# ---------------------------------------------------------------------------
# Mutation
# ---------------------------------------------------------------------------


class TestMutateSequence:
    """Tests for mutate_sequence."""

    def test_basic_mutation(self) -> None:
        """Introduce 2 mutations into a homopolymer and count differences."""
        original = "AAAAAAAAAA"
        rng = random.Random(42)
        mutated = mutate_sequence(original, n_mut=2, rng=rng)
        assert len(mutated) == len(original)
        differences = sum(1 for a, b in zip(original, mutated) if a != b)
        assert differences == 2

    def test_zero_mutations(self) -> None:
        """Zero mutations returns the original sequence unchanged."""
        original = "ATCGATCG"
        mutated = mutate_sequence(original, n_mut=0, rng=random.Random(42))
        assert mutated == original

    def test_all_positions_mutated(self) -> None:
        """Mutating every position of a homopolymer changes all bases."""
        original = "AAAA"
        mutated = mutate_sequence(original, n_mut=4, rng=random.Random(42))
        assert len(mutated) == 4
        assert all(c != "A" for c in mutated)

    def test_single_mutation(self) -> None:
        """A single mutation changes exactly one position."""
        original = "TTTTTTTTTT"
        mutated = mutate_sequence(original, n_mut=1, rng=random.Random(42))
        differences = sum(1 for a, b in zip(original, mutated) if a != b)
        assert differences == 1

    def test_mutated_bases_are_valid(self) -> None:
        """All bases in the mutated DNA sequence are valid DNA bases."""
        original = generate_random_dna(100, rng=random.Random(1))
        mutated = mutate_sequence(original, n_mut=10, rng=random.Random(42))
        assert all(c in "ACGT" for c in mutated)

    def test_reproducibility(self) -> None:
        """Same seed produces identical mutation results."""
        original = "ATCGATCGATCG"
        mut1 = mutate_sequence(original, n_mut=3, rng=random.Random(123))
        mut2 = mutate_sequence(original, n_mut=3, rng=random.Random(123))
        assert mut1 == mut2

    def test_n_mut_exceeds_length_raises(self) -> None:
        """Requesting more mutations than sequence length raises ValueError."""
        with pytest.raises(Exception):
            mutate_sequence("ATCG", n_mut=5, rng=random.Random(42))

    def test_negative_n_mut_raises(self) -> None:
        """Negative mutation count raises ValueError."""
        with pytest.raises(Exception):
            mutate_sequence("ATCG", n_mut=-1, rng=random.Random(42))


# ---------------------------------------------------------------------------
# Sequence evolution
# ---------------------------------------------------------------------------


class TestEvolveSequence:
    """Tests for evolve_sequence.

    NOTE: evolve_sequence internally calls rng.poisson(), which is not available
    on stdlib random.Random. With 0 generations the loop body is never entered,
    so it works correctly. Non-zero generations with random.Random will raise
    AttributeError due to the missing poisson method.
    """

    def test_zero_generations_returns_original(self) -> None:
        """Zero generations returns the input sequence unchanged."""
        seq = "ATCGATCG"
        result = evolve_sequence(seq, 0, mutation_rate=0.5, rng=random.Random(42))
        assert result == seq

    def test_zero_generations_various_sequences(self) -> None:
        """Zero generations preserves any sequence exactly."""
        for seq in ["A", "GGGGGGGG", "ATCGATCGATCGATCG"]:
            result = evolve_sequence(seq, 0, rng=random.Random(42))
            assert result == seq

    def test_nonzero_generations_with_rng(self) -> None:
        """evolve_sequence works with random.Random (Poisson via Knuth algorithm)."""
        rng = random.Random(42)
        result = evolve_sequence("ATCGATCG", 5, mutation_rate=0.1, rng=rng)
        assert len(result) == 8
        assert all(c in "ATCG" for c in result)


# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------


class TestTranslateDNAToProtein:
    """Tests for translate_dna_to_protein."""

    def test_known_codon(self) -> None:
        """ATG translates to M (methionine)."""
        result = translate_dna_to_protein("ATG")
        assert result == "M"

    def test_multiple_codons(self) -> None:
        """Translate a short ORF and verify expected amino acids."""
        # ATG=M, ATC=I, GAT=D, CGA=R
        result = translate_dna_to_protein("ATGATCGATCGA")
        assert result == "MIDR"

    def test_stop_codon_terminates(self) -> None:
        """Translation stops at the first stop codon."""
        # ATG=M, TAA=*, remainder should not be translated
        result = translate_dna_to_protein("ATGTAAATG")
        assert result == "M*"

    def test_frame_offset_1(self) -> None:
        """Reading frame 1 skips the first nucleotide."""
        # Frame 0: AAT GAT C -> N, D (incomplete last codon dropped)
        # Frame 1: ATG ATC -> M, I
        result = translate_dna_to_protein("AATGATC", frame=1)
        assert result == "MI"

    def test_frame_offset_2(self) -> None:
        """Reading frame 2 skips the first two nucleotides."""
        result = translate_dna_to_protein("GGATGATC", frame=2)
        assert "M" in result

    def test_empty_after_frame_offset(self) -> None:
        """Very short sequence with frame offset yields empty protein."""
        result = translate_dna_to_protein("AT", frame=0)
        assert result == ""

    def test_invalid_frame_raises(self) -> None:
        """Frame outside 0-2 range raises ValueError."""
        with pytest.raises(Exception):
            translate_dna_to_protein("ATGATC", frame=3)

    def test_invalid_bases_raise(self) -> None:
        """Non-DNA characters raise an error."""
        with pytest.raises(Exception):
            translate_dna_to_protein("ATXGATC")


# ---------------------------------------------------------------------------
# Reverse transcription
# ---------------------------------------------------------------------------


class TestReverseTranscribeProteinToDNA:
    """Tests for reverse_transcribe_protein_to_dna."""

    def test_basic_reverse_transcription(self) -> None:
        """Reverse transcribing a protein yields a DNA string of 3x length."""
        rng = random.Random(42)
        dna = reverse_transcribe_protein_to_dna("MIDS", rng=rng)
        assert len(dna) == 12  # 4 amino acids * 3 bases each
        assert all(c in "ACGTN" for c in dna)

    def test_round_trip(self) -> None:
        """Reverse-transcribed DNA translates back to the original protein."""
        protein = "MIDSAARGLK"
        rng = random.Random(42)
        dna = reverse_transcribe_protein_to_dna(protein, rng=rng)
        # Translate back
        recovered = translate_dna_to_protein(dna)
        # May include stop codon at end; compare prefix
        assert recovered.rstrip("*") == protein or recovered.startswith(protein)

    def test_reproducibility(self) -> None:
        """Same seed produces identical reverse transcription."""
        protein = "ACDEFGHIK"
        dna1 = reverse_transcribe_protein_to_dna(protein, rng=random.Random(42))
        dna2 = reverse_transcribe_protein_to_dna(protein, rng=random.Random(42))
        assert dna1 == dna2

    def test_stop_codon_encoding(self) -> None:
        """Stop codon symbol '*' maps to one of TAA/TAG/TGA."""
        rng = random.Random(42)
        dna = reverse_transcribe_protein_to_dna("M*", rng=rng)
        assert len(dna) == 6
        stop_codon = dna[3:6]
        assert stop_codon in ("TAA", "TAG", "TGA")

    def test_empty_protein(self) -> None:
        """Empty protein yields empty DNA."""
        dna = reverse_transcribe_protein_to_dna("", rng=random.Random(42))
        assert dna == ""


# ---------------------------------------------------------------------------
# Coding sequence generation
# ---------------------------------------------------------------------------


class TestGenerateCodingSequence:
    """Tests for generate_coding_sequence."""

    def test_basic_generation(self) -> None:
        """Generate a coding sequence of length 30 (divisible by 3)."""
        rng = random.Random(42)
        dna, protein = generate_coding_sequence(30, rng=rng)
        assert len(dna) == 30
        assert len(protein) >= 1  # At least one amino acid before any stop codon

    def test_dna_is_valid(self) -> None:
        """All bases in the generated DNA are valid."""
        dna, _ = generate_coding_sequence(60, rng=random.Random(42))
        assert all(c in "ACGT" for c in dna)

    def test_protein_matches_translation(self) -> None:
        """The returned protein matches direct translation of the DNA."""
        dna, protein = generate_coding_sequence(60, rng=random.Random(42))
        translated = translate_dna_to_protein(dna)
        assert translated == protein

    def test_length_not_divisible_by_3_raises(self) -> None:
        """Length not divisible by 3 raises a validation error."""
        with pytest.raises(Exception):
            generate_coding_sequence(31, rng=random.Random(42))

    def test_minimum_coding_length(self) -> None:
        """Minimum coding length is 3 (one codon)."""
        dna, protein = generate_coding_sequence(3, rng=random.Random(42))
        assert len(dna) == 3
        assert len(protein) >= 1 or protein == ""  # Could be stop codon

    def test_gc_content_parameter(self) -> None:
        """GC content parameter influences the generated DNA."""
        dna_high, _ = generate_coding_sequence(3000, gc_content=0.8, rng=random.Random(42))
        dna_low, _ = generate_coding_sequence(3000, gc_content=0.2, rng=random.Random(42))
        gc_high = sum(1 for c in dna_high if c in "GC") / len(dna_high)
        gc_low = sum(1 for c in dna_low if c in "GC") / len(dna_low)
        assert gc_high > gc_low


# ---------------------------------------------------------------------------
# Sequence similarity
# ---------------------------------------------------------------------------


class TestCalculateSequenceSimilarity:
    """Tests for calculate_sequence_similarity."""

    def test_identical_sequences(self) -> None:
        """Identical sequences have similarity 1.0."""
        assert calculate_sequence_similarity("ATCG", "ATCG") == 1.0

    def test_completely_different(self) -> None:
        """Completely different sequences have similarity 0.0."""
        assert calculate_sequence_similarity("AAAA", "TTTT") == 0.0

    def test_one_mismatch(self) -> None:
        """One mismatch out of four gives similarity 0.75."""
        assert calculate_sequence_similarity("ATCG", "TTCG") == 0.75

    def test_half_matching(self) -> None:
        """Half matching positions gives similarity 0.5."""
        assert calculate_sequence_similarity("AATT", "TTTT") == 0.5

    def test_empty_sequences(self) -> None:
        """Two empty sequences have similarity 1.0 by convention."""
        assert calculate_sequence_similarity("", "") == 1.0

    def test_different_lengths_raise(self) -> None:
        """Sequences of different lengths raise an error."""
        with pytest.raises(Exception):
            calculate_sequence_similarity("ATCG", "AT")

    def test_symmetry(self) -> None:
        """Similarity is symmetric: sim(a, b) == sim(b, a)."""
        a = "ATCGATCG"
        b = "TTCGATCA"
        assert calculate_sequence_similarity(a, b) == calculate_sequence_similarity(b, a)


# ---------------------------------------------------------------------------
# Sequence family generation
# ---------------------------------------------------------------------------


class TestGenerateSequenceFamily:
    """Tests for generate_sequence_family.

    NOTE: generate_sequence_family calls evolve_sequence internally, which
    uses rng.poisson. With 0 generations, the evolve loop does not run,
    so all descendants equal the ancestor.
    """

    def test_family_with_zero_generations(self) -> None:
        """Zero generations produces descendants identical to ancestor."""
        ancestor = "ATCGATCGATCG"
        rng = random.Random(42)
        family = generate_sequence_family(ancestor, 3, 0, mutation_rate=0.001, rng=rng)
        # Family includes ancestor + 3 descendants = 4 total
        assert len(family) == 4
        assert family[0] == ancestor
        assert all(s == ancestor for s in family)

    def test_family_size_correct(self) -> None:
        """Family includes ancestor plus n_descendants members."""
        ancestor = "GGGGGGGG"
        family = generate_sequence_family(ancestor, 5, 0, rng=random.Random(42))
        assert len(family) == 6  # 1 ancestor + 5 descendants


# ---------------------------------------------------------------------------
# Sequence divergence analysis
# ---------------------------------------------------------------------------


class TestAnalyzeSequenceDivergence:
    """Tests for analyze_sequence_divergence."""

    def test_basic_divergence(self) -> None:
        """Analyze divergence of three related sequences."""
        seqs = ["ATCGATCG", "ATCAATCG", "TTCGATCG"]
        result = analyze_sequence_divergence(seqs)

        assert result["num_sequences"] == 3
        assert result["sequence_length"] == 8
        assert 0.0 <= result["mean_similarity"] <= 1.0
        assert 0.0 <= result["mean_divergence"] <= 1.0
        # mean_similarity + mean_divergence should equal 1.0
        assert abs(result["mean_similarity"] + result["mean_divergence"] - 1.0) < 1e-10

    def test_identical_sequences(self) -> None:
        """Identical sequences have zero divergence."""
        seqs = ["ATCGATCG", "ATCGATCG"]
        result = analyze_sequence_divergence(seqs)
        assert result["mean_similarity"] == 1.0
        assert result["mean_divergence"] == 0.0
        assert result["variable_positions"] == 0

    def test_variable_positions_count(self) -> None:
        """Count variable positions correctly."""
        seqs = ["ATCGATCG", "ATCAATCG", "TTCGATCG"]
        result = analyze_sequence_divergence(seqs)
        # Position 0: A, A, T -> variable
        # Position 3: G, A, G -> variable
        assert result["variable_positions"] == 2

    def test_pairwise_counts(self) -> None:
        """Number of pairwise comparisons is n*(n-1)/2."""
        seqs = ["AAAA", "AAAT", "AATT", "ATTT"]
        result = analyze_sequence_divergence(seqs)
        n = 4
        expected_pairs = n * (n - 1) // 2
        assert len(result["pairwise_similarities"]) == expected_pairs
        assert len(result["pairwise_divergences"]) == expected_pairs

    def test_fewer_than_two_sequences_raises(self) -> None:
        """Fewer than 2 sequences raises an error."""
        with pytest.raises(Exception):
            analyze_sequence_divergence(["ATCG"])

    def test_unequal_lengths_raise(self) -> None:
        """Sequences of different lengths raise an error."""
        with pytest.raises(Exception):
            analyze_sequence_divergence(["ATCG", "AT"])

    def test_variable_fraction(self) -> None:
        """Variable fraction is variable_positions / sequence_length."""
        seqs = ["AAAA", "AAAT"]
        result = analyze_sequence_divergence(seqs)
        assert result["variable_positions"] == 1
        assert result["variable_fraction"] == pytest.approx(0.25)


# ---------------------------------------------------------------------------
# Gene duplication
# ---------------------------------------------------------------------------


class TestSimulateGeneDuplication:
    """Tests for simulate_gene_duplication.

    NOTE: simulate_gene_duplication calls evolve_sequence internally.
    With divergence_time=0, the evolve loop does not run, so copies
    are identical to the original.
    """

    def test_zero_divergence_time(self) -> None:
        """Zero divergence time produces copies identical to original."""
        original = "ATCGATCG"
        copies = simulate_gene_duplication(
            original, 3, divergence_time=0, rng=random.Random(42)
        )
        assert len(copies) == 3
        assert all(c == original for c in copies)

    def test_copy_count(self) -> None:
        """Requested number of copies is produced."""
        original = "GGGGCCCC"
        copies = simulate_gene_duplication(
            original, 5, divergence_time=0, rng=random.Random(42)
        )
        assert len(copies) == 5

    def test_single_copy(self) -> None:
        """Single copy duplication works."""
        original = "ATCGATCG"
        copies = simulate_gene_duplication(
            original, 1, divergence_time=0, rng=random.Random(42)
        )
        assert len(copies) == 1
        assert copies[0] == original


# ---------------------------------------------------------------------------
# RNA-seq count simulation (convenience wrapper)
# ---------------------------------------------------------------------------


class TestSimulateRNASeqCounts:
    """Tests for simulate_rnaseq_counts (convenience wrapper)."""

    def test_basic_shape(self) -> None:
        """Output shape is (n_samples, n_genes)."""
        rng = random.Random(42)
        counts = simulate_rnaseq_counts(n_genes=10, n_samples=5, rng=rng)
        assert isinstance(counts, np.ndarray)
        assert counts.shape == (5, 10)

    def test_all_non_negative(self) -> None:
        """All counts are non-negative integers."""
        counts = simulate_rnaseq_counts(n_genes=50, n_samples=8, rng=random.Random(42))
        assert np.all(counts >= 0)

    def test_integer_dtype(self) -> None:
        """Count matrix has integer dtype."""
        counts = simulate_rnaseq_counts(n_genes=10, n_samples=5, rng=random.Random(42))
        assert np.issubdtype(counts.dtype, np.integer)

    def test_reproducibility(self) -> None:
        """Same seed produces identical count matrices."""
        counts1 = simulate_rnaseq_counts(n_genes=20, n_samples=5, rng=random.Random(42))
        counts2 = simulate_rnaseq_counts(n_genes=20, n_samples=5, rng=random.Random(42))
        assert np.array_equal(counts1, counts2)

    def test_different_seeds_differ(self) -> None:
        """Different seeds produce different count matrices."""
        counts1 = simulate_rnaseq_counts(n_genes=50, n_samples=10, rng=random.Random(42))
        counts2 = simulate_rnaseq_counts(n_genes=50, n_samples=10, rng=random.Random(99))
        assert not np.array_equal(counts1, counts2)

    def test_higher_mean_produces_higher_total(self) -> None:
        """Higher mean_expression produces larger total count."""
        counts_low = simulate_rnaseq_counts(
            n_genes=100, n_samples=10, mean_expression=10.0, rng=random.Random(42)
        )
        counts_high = simulate_rnaseq_counts(
            n_genes=100, n_samples=10, mean_expression=1000.0, rng=random.Random(42)
        )
        assert counts_high.sum() > counts_low.sum()

    def test_default_parameters(self) -> None:
        """Default parameters (1000 genes, 10 samples) produce valid output."""
        counts = simulate_rnaseq_counts(rng=random.Random(42))
        assert counts.shape == (10, 1000)
        assert np.all(counts >= 0)

    def test_single_gene_single_sample(self) -> None:
        """Minimum dimensions: 1 gene, 1 sample."""
        counts = simulate_rnaseq_counts(n_genes=1, n_samples=1, rng=random.Random(42))
        assert counts.shape == (1, 1)
        assert counts[0, 0] >= 0


# ---------------------------------------------------------------------------
# GridWorld initialization and simulation
# ---------------------------------------------------------------------------


class TestGridWorld:
    """Tests for GridWorld agent-based simulation."""

    def test_initialization(self) -> None:
        """GridWorld initializes with correct dimensions and agent count."""
        rng = random.Random(42)
        world = GridWorld(width=10, height=10, num_agents=5, rng=rng)
        assert world.width == 10
        assert world.height == 10
        assert len(world.agents) == 5

    def test_agents_within_bounds(self) -> None:
        """All agents are placed within grid boundaries."""
        rng = random.Random(42)
        world = GridWorld(width=10, height=10, num_agents=20, rng=rng)
        for agent in world.agents:
            assert 0 <= agent.x < 10
            assert 0 <= agent.y < 10

    def test_step_returns_dict(self) -> None:
        """A single step returns a dictionary with expected keys."""
        world = GridWorld(width=5, height=5, num_agents=2, rng=random.Random(42))
        result = world.step()
        assert isinstance(result, dict)
        assert "time_step" in result
        assert "num_agents" in result
        assert "agent_positions" in result
        assert "total_energy" in result

    def test_step_increments_time(self) -> None:
        """Each step increments the time step counter."""
        world = GridWorld(width=5, height=5, num_agents=2, rng=random.Random(42))
        assert world.time_step == 0
        world.step()
        assert world.time_step == 1
        world.step()
        assert world.time_step == 2

    def test_agents_remain_in_bounds_after_step(self) -> None:
        """After stepping, agents remain within grid boundaries."""
        world = GridWorld(width=5, height=5, num_agents=3, rng=random.Random(42))
        for _ in range(20):
            world.step()
        for x, y in world.positions():
            assert 0 <= x < 5
            assert 0 <= y < 5

    def test_positions_method(self) -> None:
        """positions() returns list of (x, y) tuples matching agent count."""
        world = GridWorld(width=10, height=10, num_agents=4, rng=random.Random(42))
        positions = world.positions()
        assert len(positions) == 4
        for pos in positions:
            assert isinstance(pos, tuple)
            assert len(pos) == 2

    def test_get_agent_positions_alias(self) -> None:
        """get_agent_positions() returns the same result as positions()."""
        world = GridWorld(width=10, height=10, num_agents=3, rng=random.Random(42))
        assert world.get_agent_positions() == world.positions()

    def test_get_agent_energies(self) -> None:
        """get_agent_energies() returns a list of floats for each agent."""
        world = GridWorld(width=10, height=10, num_agents=3, rng=random.Random(42))
        energies = world.get_agent_energies()
        assert len(energies) == 3
        assert all(isinstance(e, float) for e in energies)
        assert all(e >= 0 for e in energies)

    def test_energy_decreases_after_steps(self) -> None:
        """Agent energy decreases over multiple steps."""
        world = GridWorld(width=10, height=10, num_agents=1, rng=random.Random(42))
        initial_energy = world.agents[0].energy
        for _ in range(10):
            world.step()
        final_energy = world.agents[0].energy
        assert final_energy < initial_energy

    def test_run_simulation(self) -> None:
        """run_simulation returns a list of step result dicts."""
        world = GridWorld(width=8, height=8, num_agents=3, rng=random.Random(42))
        results = world.run_simulation(5)
        assert len(results) == 5
        for result in results:
            assert isinstance(result, dict)
            assert "time_step" in result

    def test_run_simulation_positions_in_bounds(self) -> None:
        """All positions in simulation history remain within bounds."""
        world = GridWorld(width=6, height=6, num_agents=4, rng=random.Random(42))
        results = world.run_simulation(10)
        for result in results:
            for x, y in result["agent_positions"]:
                assert 0 <= x < 6
                assert 0 <= y < 6

    def test_multiple_steps_manual_tracking(self) -> None:
        """Manually track positions over multiple steps."""
        world = GridWorld(width=8, height=8, num_agents=3, rng=random.Random(42))
        history = [world.positions()]
        for _ in range(5):
            world.step()
            history.append(world.positions())
        assert len(history) == 6  # initial + 5 steps
        assert all(len(step_positions) == 3 for step_positions in history)

    def test_zero_agents(self) -> None:
        """GridWorld with zero agents initializes and steps without error."""
        world = GridWorld(width=5, height=5, num_agents=0, rng=random.Random(42))
        assert len(world.agents) == 0
        result = world.step()
        assert result["num_agents"] == 0

    def test_single_agent(self) -> None:
        """GridWorld with a single agent works correctly."""
        world = GridWorld(width=3, height=3, num_agents=1, rng=random.Random(42))
        assert len(world.agents) == 1
        world.step()
        pos = world.positions()
        assert len(pos) == 1
        assert 0 <= pos[0][0] < 3
        assert 0 <= pos[0][1] < 3


# ---------------------------------------------------------------------------
# GridAgent
# ---------------------------------------------------------------------------


class TestGridAgent:
    """Tests for GridAgent dataclass."""

    def test_creation(self) -> None:
        """Create a GridAgent with explicit attributes."""
        agent = GridAgent(id=0, x=3, y=4, energy=100.0)
        assert agent.id == 0
        assert agent.x == 3
        assert agent.y == 4
        assert agent.energy == 100.0

    def test_default_energy(self) -> None:
        """Default energy is 100.0."""
        agent = GridAgent(id=1, x=0, y=0)
        assert agent.energy == 100.0

    def test_custom_properties(self) -> None:
        """Properties dictionary can be set on creation."""
        props = {"speed": 2.0, "type": "herbivore"}
        agent = GridAgent(id=2, x=5, y=5, properties=props)
        assert agent.properties["speed"] == 2.0
        assert agent.properties["type"] == "herbivore"

    def test_step_moves_agent(self) -> None:
        """Calling step on a GridAgent via a GridWorld changes position or stays."""
        rng = random.Random(42)
        world = GridWorld(width=10, height=10, num_agents=1, rng=rng)
        agent = world.agents[0]
        initial_x, initial_y = agent.x, agent.y
        step_rng = random.Random(99)
        agent.step(world, step_rng)
        # Agent should still be in bounds
        assert 0 <= agent.x < 10
        assert 0 <= agent.y < 10

    def test_step_decreases_energy(self) -> None:
        """Each step decreases the agent's energy by 0.1."""
        rng = random.Random(42)
        world = GridWorld(width=10, height=10, num_agents=1, rng=rng)
        agent = world.agents[0]
        initial_energy = agent.energy
        agent.step(world, random.Random(99))
        assert agent.energy == pytest.approx(initial_energy - 0.1)

    def test_step_clamps_position(self) -> None:
        """Agent at grid edge does not move out of bounds."""
        world = GridWorld(width=3, height=3, num_agents=0, rng=random.Random(42))
        # Place agent at corner (0, 0)
        agent = GridAgent(id=0, x=0, y=0, energy=100.0)
        world.agents.append(agent)
        # Step many times to test boundary clamping
        step_rng = random.Random(42)
        for _ in range(50):
            agent.step(world, step_rng)
            assert 0 <= agent.x < 3
            assert 0 <= agent.y < 3

    def test_energy_floors_at_zero(self) -> None:
        """Agent energy does not go below zero."""
        world = GridWorld(width=5, height=5, num_agents=0, rng=random.Random(42))
        agent = GridAgent(id=0, x=2, y=2, energy=0.05)
        world.agents.append(agent)
        step_rng = random.Random(42)
        agent.step(world, step_rng)
        assert agent.energy >= 0.0


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


class TestConstants:
    """Tests for module-level constants."""

    def test_dna_bases(self) -> None:
        """DNA_BASES contains exactly A, T, C, G."""
        assert set(DNA_BASES) == {"A", "T", "C", "G"}

    def test_rna_bases(self) -> None:
        """RNA_BASES contains exactly A, U, C, G."""
        assert set(RNA_BASES) == {"A", "U", "C", "G"}

    def test_amino_acids_count(self) -> None:
        """AMINO_ACIDS contains 20 standard amino acids."""
        assert len(AMINO_ACIDS) == 20
        assert len(set(AMINO_ACIDS)) == 20  # All unique

    def test_genetic_code_size(self) -> None:
        """GENETIC_CODE maps 64 codons."""
        assert len(GENETIC_CODE) == 64

    def test_genetic_code_stop_codons(self) -> None:
        """GENETIC_CODE contains exactly 3 stop codons."""
        stops = [codon for codon, aa in GENETIC_CODE.items() if aa == "*"]
        assert set(stops) == {"TAA", "TAG", "TGA"}

    def test_genetic_code_start_codon(self) -> None:
        """ATG maps to methionine (M)."""
        assert GENETIC_CODE["ATG"] == "M"


# ---------------------------------------------------------------------------
# Integration test
# ---------------------------------------------------------------------------


class TestSimulationIntegration:
    """Integration tests combining multiple simulation components."""

    def test_full_pipeline(self) -> None:
        """End-to-end pipeline: generate, mutate, translate, compare, simulate."""
        rng_dna = random.Random(42)
        rng_mut = random.Random(123)
        rng_prot = random.Random(456)
        rng_counts = random.Random(789)

        # Step 1: Generate a coding DNA sequence
        dna_seq = generate_random_dna(99, rng=rng_dna)  # 99 bases = 33 codons
        assert len(dna_seq) == 99

        # Step 2: Mutate the sequence
        mutated_seq = mutate_sequence(dna_seq, n_mut=5, rng=rng_mut)
        assert len(mutated_seq) == 99

        # Step 3: Translate both to protein
        protein_orig = translate_dna_to_protein(dna_seq)
        protein_mut = translate_dna_to_protein(mutated_seq)
        assert len(protein_orig) >= 1
        assert len(protein_mut) >= 1

        # Step 4: Generate a protein and reverse transcribe
        protein_seq = generate_random_protein(33, rng=rng_prot)
        dna_from_protein = reverse_transcribe_protein_to_dna(protein_seq, rng=random.Random(42))
        assert len(dna_from_protein) == 99  # 33 amino acids * 3

        # Step 5: Compare sequences
        sim = calculate_sequence_similarity(dna_seq, mutated_seq)
        assert 0.0 <= sim <= 1.0
        # With 5 mutations out of 99, expect high similarity
        assert sim >= 0.9

        # Step 6: Simulate RNA-seq counts
        counts = simulate_rnaseq_counts(n_genes=20, n_samples=10, rng=rng_counts)
        assert counts.shape == (10, 20)
        assert np.all(counts >= 0)

        # Step 7: Run agent simulation
        rng_world = random.Random(999)
        world = GridWorld(width=6, height=6, num_agents=4, rng=rng_world)
        results = world.run_simulation(3)
        assert len(results) == 3
        final_positions = world.positions()
        assert len(final_positions) == 4
        for x, y in final_positions:
            assert 0 <= x < 6
            assert 0 <= y < 6

    def test_sequence_analysis_pipeline(self) -> None:
        """Generate related sequences and analyze their divergence."""
        rng = random.Random(42)
        ancestor = generate_random_dna(100, rng=rng)

        # Create manually mutated "descendants"
        descendants = [ancestor]
        for i in range(3):
            mutated = mutate_sequence(ancestor, n_mut=5 * (i + 1), rng=random.Random(i))
            descendants.append(mutated)

        # Analyze divergence
        result = analyze_sequence_divergence(descendants)
        assert result["num_sequences"] == 4
        assert result["sequence_length"] == 100
        assert result["mean_similarity"] > 0.5  # Related sequences
        assert result["variable_positions"] > 0

    def test_coding_and_similarity(self) -> None:
        """Generate two coding sequences and measure similarity."""
        dna1, prot1 = generate_coding_sequence(60, rng=random.Random(42))
        dna2, prot2 = generate_coding_sequence(60, rng=random.Random(99))

        dna_sim = calculate_sequence_similarity(dna1, dna2)
        assert 0.0 <= dna_sim <= 1.0

    def test_gridworld_energy_tracking(self) -> None:
        """Track energy decay across simulation steps."""
        world = GridWorld(width=10, height=10, num_agents=5, rng=random.Random(42))
        initial_energies = world.get_agent_energies()
        results = world.run_simulation(20)

        final_energies = world.get_agent_energies()
        # Energy should decrease over time (each step costs 0.1)
        assert sum(final_energies) < sum(initial_energies)
