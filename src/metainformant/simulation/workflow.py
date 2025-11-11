"""Simulation workflow orchestration for METAINFORMANT.

Provides end-to-end simulation pipelines for synthetic data generation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..core import io, paths, validation
from ..core.logging import get_logger

logger = get_logger(__name__)


def run_sequence_simulation_workflow(
    output_dir: Path | str,
    *,
    n_sequences: int = 100,
    sequence_length: int = 1000,
    gc_content: float = 0.5,
    seed: int | None = None,
) -> dict[str, Any]:
    """Run end-to-end sequence simulation workflow.
    
    Args:
        output_dir: Directory for output files
        n_sequences: Number of sequences to generate
        sequence_length: Length of each sequence
        gc_content: GC content for sequences
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with workflow results and metadata
    """
    from . import sequences
    import random
    
    validation.validate_type(n_sequences, int, "n_sequences")
    validation.validate_range(n_sequences, min_val=1, name="n_sequences")
    validation.validate_range(sequence_length, min_val=1, name="sequence_length")
    validation.validate_range(gc_content, min_val=0.0, max_val=1.0, name="gc_content")
    
    output_dir = Path(output_dir)
    paths.ensure_directory(output_dir)
    
    logger.info(f"Starting sequence simulation workflow: {n_sequences} sequences")
    
    # Setup RNG
    rng = random.Random(seed) if seed is not None else random.Random()
    
    # Generate sequences
    sequences_list = []
    for i in range(n_sequences):
        seq = sequences.generate_random_dna(sequence_length, gc_content=gc_content, rng=rng)
        sequences_list.append(seq)
    
    # Save results
    results = {
        "n_sequences": n_sequences,
        "sequence_length": sequence_length,
        "gc_content": gc_content,
        "seed": seed,
        "sequences": sequences_list
    }
    
    output_file = output_dir / "sequences.json"
    io.dump_json(results, output_file)
    
    logger.info(f"Sequence simulation complete. Results saved to {output_file}")
    
    return {
        "status": "completed",
        "output_file": str(output_file),
        "n_sequences": n_sequences
    }


def run_agent_simulation_workflow(
    output_dir: Path | str,
    *,
    width: int = 50,
    height: int = 50,
    num_agents: int = 100,
    num_steps: int = 100,
    seed: int | None = None,
) -> dict[str, Any]:
    """Run end-to-end agent-based simulation workflow.
    
    Args:
        output_dir: Directory for output files
        width: Grid width
        height: Grid height
        num_agents: Number of agents
        num_steps: Number of simulation steps
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with workflow results and metadata
    """
    from . import agents
    import random
    
    validation.validate_range(width, min_val=1, name="width")
    validation.validate_range(height, min_val=1, name="height")
    validation.validate_range(num_agents, min_val=0, name="num_agents")
    validation.validate_range(num_steps, min_val=0, name="num_steps")
    
    output_dir = Path(output_dir)
    paths.ensure_directory(output_dir)
    
    logger.info(f"Starting agent simulation workflow: {num_agents} agents, {num_steps} steps")
    
    # Setup RNG
    rng = random.Random(seed) if seed is not None else random.Random()
    
    # Create world
    world = agents.GridWorld(width, height, num_agents, rng=rng)
    
    # Run simulation
    position_history = []
    for step in range(num_steps):
        world.step()
        if step % 10 == 0:  # Record every 10 steps
            position_history.append({
                "step": step,
                "positions": world.positions()
            })
    
    # Save results
    results = {
        "width": width,
        "height": height,
        "num_agents": num_agents,
        "num_steps": num_steps,
        "seed": seed,
        "final_positions": world.positions(),
        "position_history": position_history
    }
    
    output_file = output_dir / "agent_simulation.json"
    io.dump_json(results, output_file)
    
    logger.info(f"Agent simulation complete. Results saved to {output_file}")
    
    return {
        "status": "completed",
        "output_file": str(output_file),
        "num_steps": num_steps,
        "num_agents": num_agents
    }


def run_popgen_simulation_workflow(
    output_dir: Path | str,
    *,
    n_sequences: int = 50,
    sequence_length: int = 1000,
    nucleotide_diversity: float = 0.01,
    seed: int | None = None,
) -> dict[str, Any]:
    """Run end-to-end population genetics simulation workflow.
    
    Args:
        output_dir: Directory for output files
        n_sequences: Number of sequences
        sequence_length: Length of each sequence
        nucleotide_diversity: Target nucleotide diversity
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with workflow results and metadata
    """
    from . import popgen
    import random
    
    validation.validate_range(n_sequences, min_val=1, name="n_sequences")
    validation.validate_range(sequence_length, min_val=1, name="sequence_length")
    validation.validate_range(nucleotide_diversity, min_val=0.0, name="nucleotide_diversity")
    
    output_dir = Path(output_dir)
    paths.ensure_directory(output_dir)
    
    logger.info(f"Starting popgen simulation workflow: {n_sequences} sequences")
    
    # Setup RNG
    rng = random.Random(seed) if seed is not None else random.Random()
    
    # Generate population
    sequences_list = popgen.generate_population_sequences(
        n_sequences,
        sequence_length,
        nucleotide_diversity=nucleotide_diversity,
        rng=rng
    )
    
    # Save results
    results = {
        "n_sequences": n_sequences,
        "sequence_length": sequence_length,
        "nucleotide_diversity": nucleotide_diversity,
        "seed": seed,
        "sequences": sequences_list
    }
    
    output_file = output_dir / "popgen_simulation.json"
    io.dump_json(results, output_file)
    
    logger.info(f"Popgen simulation complete. Results saved to {output_file}")
    
    return {
        "status": "completed",
        "output_file": str(output_file),
        "n_sequences": n_sequences
    }


