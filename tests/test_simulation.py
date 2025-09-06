from __future__ import annotations

import random

from metainformant.simulation import (
    Agent,
    GridWorld,
    generate_random_dna,
    generate_random_protein,
    mutate_sequence,
    simulate_counts_negative_binomial,
)


def test_generate_random_dna_basic():
    """Test DNA sequence generation with proper length and nucleotides."""
    rng = random.Random(42)
    seq = generate_random_dna(100, rng=rng)
    assert len(seq) == 100
    assert all(c in "ACGT" for c in seq)

    # Test reproducibility
    rng2 = random.Random(42)
    seq2 = generate_random_dna(100, rng=rng2)
    assert seq == seq2


def test_generate_random_dna_edge_cases():
    """Test DNA generation edge cases."""
    rng = random.Random(42)
    # Empty sequence
    seq = generate_random_dna(0, rng=rng)
    assert len(seq) == 0
    assert seq == ""

    # Single nucleotide
    seq = generate_random_dna(1, rng=rng)
    assert len(seq) == 1
    assert seq in "ACGT"


def test_mutate_sequence_basic():
    """Test sequence mutation functionality."""
    original = "AAAAAAAAAA"  # 10 A's
    rng = random.Random(42)
    mutated = mutate_sequence(original, n_mut=2, rng=rng)
    assert len(mutated) == len(original)

    # Exactly 2 mutations should occur
    differences = sum(1 for a, b in zip(original, mutated) if a != b)
    assert differences == 2


def test_mutate_sequence_reproducible():
    """Test mutation reproducibility with seeds."""
    original = "ATCGATCGATCG"
    rng1 = random.Random(123)
    rng2 = random.Random(123)
    mut1 = mutate_sequence(original, n_mut=3, rng=rng1)
    mut2 = mutate_sequence(original, n_mut=3, rng=rng2)
    assert mut1 == mut2


def test_mutate_sequence_edge_cases():
    """Test mutation edge cases."""
    # No mutations
    original = "ATCG"
    rng = random.Random(42)
    mutated = mutate_sequence(original, n_mut=0, rng=rng)
    assert mutated == original

    # More mutations than sequence length
    mutated = mutate_sequence(original, n_mut=10, rng=rng)
    assert len(mutated) == len(original)


def test_generate_random_protein_basic():
    """Test protein sequence generation."""
    rng = random.Random(42)
    seq = generate_random_protein(50, rng=rng)
    assert len(seq) == 50
    # Check that all characters are valid amino acids (20 standard)
    valid_aa = "ACDEFGHIKLMNPQRSTVWY"
    assert all(c in valid_aa for c in seq)

    # Test reproducibility
    rng2 = random.Random(42)
    seq2 = generate_random_protein(50, rng=rng2)
    assert seq == seq2


def test_generate_random_protein_edge_cases():
    """Test protein generation edge cases."""
    rng = random.Random(42)
    seq = generate_random_protein(0, rng=rng)
    assert len(seq) == 0

    seq = generate_random_protein(1, rng=rng)
    assert len(seq) == 1


def test_simulate_counts_negative_binomial():
    """Test negative binomial count simulation."""
    rng = random.Random(42)
    counts = simulate_counts_negative_binomial(num_genes=10, num_samples=5, rng=rng)
    assert len(counts) == 10  # 10 genes (rows)
    assert all(len(row) == 5 for row in counts)  # 5 samples per gene

    # All counts should be non-negative integers
    for row in counts:
        for count in row:
            assert isinstance(count, int)
            assert count >= 0


def test_simulate_counts_reproducible():
    """Test count simulation reproducibility."""
    rng1 = random.Random(123)
    rng2 = random.Random(123)
    counts1 = simulate_counts_negative_binomial(num_genes=5, num_samples=3, rng=rng1)
    counts2 = simulate_counts_negative_binomial(num_genes=5, num_samples=3, rng=rng2)
    assert counts1 == counts2


def test_simulate_counts_parameters():
    """Test count simulation with different parameters."""
    rng = random.Random(42)
    # Higher mean should generally produce higher counts
    low_counts = simulate_counts_negative_binomial(num_genes=5, num_samples=5, mean_expression=10.0, rng=rng)
    rng = random.Random(42)  # Reset for fair comparison
    high_counts = simulate_counts_negative_binomial(num_genes=5, num_samples=5, mean_expression=1000.0, rng=rng)

    low_total = sum(sum(row) for row in low_counts)
    high_total = sum(sum(row) for row in high_counts)
    assert high_total > low_total


def test_gridworld_initialization():
    """Test GridWorld agent-based simulation initialization."""
    rng = random.Random(42)
    world = GridWorld(width=10, height=10, num_agents=5, rng=rng)
    assert world.width == 10
    assert world.height == 10
    assert len(world.agents) == 5

    # Ensure agents are placed within bounds
    for agent in world.agents:
        assert 0 <= agent.x < 10
        assert 0 <= agent.y < 10


def test_gridworld_step():
    """Test GridWorld simulation step functionality."""
    rng = random.Random(42)
    world = GridWorld(width=5, height=5, num_agents=2, rng=rng)
    initial_positions = world.positions()

    # Run a simulation step
    world.step()
    new_positions = world.positions()

    # Agents should still be within bounds
    for x, y in new_positions:
        assert 0 <= x < 5
        assert 0 <= y < 5

    # Check that positions are recorded correctly
    assert len(new_positions) == 2


def test_gridworld_multiple_steps():
    """Test running multiple GridWorld simulation steps."""
    rng = random.Random(42)
    world = GridWorld(width=8, height=8, num_agents=3, rng=rng)

    # Record history manually by running multiple steps
    history = [world.positions()]  # Initial positions
    for _ in range(5):
        world.step()
        history.append(world.positions())

    assert len(history) == 6  # Initial state + 5 steps
    assert len(history[0]) == 3  # 3 agents

    # Check that each step records agent positions within bounds
    for step in history:
        for x, y in step:
            assert 0 <= x < 8
            assert 0 <= y < 8


def test_agent_creation():
    """Test individual Agent creation."""
    agent = Agent(x=3, y=4)
    assert agent.x == 3
    assert agent.y == 4


def test_agent_step():
    """Test agent movement during step."""
    rng = random.Random(42)
    world = GridWorld(width=10, height=10, num_agents=1, rng=rng)
    agent = world.agents[0]
    initial_x, initial_y = agent.x, agent.y

    # Step the agent (using world's step method)
    agent.step(world, rng)

    # Agent should still be within bounds (with wrapping)
    assert 0 <= agent.x < 10
    assert 0 <= agent.y < 10


def test_simulation_integration():
    """Integration test combining multiple simulation components."""
    rng1 = random.Random(42)
    rng2 = random.Random(123)
    rng3 = random.Random(456)
    rng4 = random.Random(789)

    # Generate synthetic data
    dna_seq = generate_random_dna(100, rng=rng1)
    mutated_seq = mutate_sequence(dna_seq, n_mut=5, rng=rng2)
    protein_seq = generate_random_protein(33, rng=rng3)  # ~100/3
    counts = simulate_counts_negative_binomial(num_genes=20, num_samples=10, rng=rng4)

    # Basic validation that all components work together
    assert len(dna_seq) == 100
    assert len(mutated_seq) == 100
    assert len(protein_seq) == 33
    assert len(counts) == 20
    assert all(len(row) == 10 for row in counts)

    # Run a small agent simulation
    rng5 = random.Random(999)
    world = GridWorld(width=6, height=6, num_agents=4, rng=rng5)

    # Manually run simulation steps
    for _ in range(3):
        world.step()

    # Verify final state
    final_positions = world.positions()
    assert len(final_positions) == 4
    for x, y in final_positions:
        assert 0 <= x < 6
        assert 0 <= y < 6
