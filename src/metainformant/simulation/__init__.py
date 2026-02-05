"""Evolutionary and Population Genetics Simulation module for METAINFORMANT.

This module provides comprehensive tools for simulating evolutionary processes,
including agent-based modeling, population genetics simulation, molecular evolution,
and sequence generation.

Capabilities:
    - **Sequence Simulation**: Generate random DNA/protein sequences and apply mutations
    - **RNA-seq Simulation**: Simulate count data using negative binomial distributions
    - **Agent-based Models**: Create agents in grid-based environments for ecological modeling
    - **Population Genetics**: Forward and backward-time coalescent simulations
    - **Molecular Evolution**: Simulate substitution models and phylogenetic trees
    - **Workflow Orchestration**: Configure and run simulation pipelines

Key Classes:
    - Agent: Individual agent for ecosystem-based modeling (typed: producer/consumer/decomposer)
    - GridAgent: Lightweight agent for grid-based spatial simulations
    - GridWorld: 2D grid environment for agent simulations
    - Ecosystem: Complex ecosystem with typed agents
    - SimulationConfig: Configuration for simulation workflows

Key Functions:
    - generate_random_dna: Generate random DNA sequences
    - generate_random_protein: Generate random protein sequences
    - mutate_sequence: Apply mutations to sequences
    - evolve_sequence: Simulate sequence evolution over generations
    - simulate_rnaseq_counts: Simulate RNA-seq count data (convenience wrapper)
    - simulate_counts_negative_binomial: Core NB count simulation
    - simulate_bulk_rnaseq: Bulk RNA-seq simulation with library sizes
    - simulate_single_cell_rnaseq: Single-cell RNA-seq with dropout
    - generate_population_sequences: Generate population of sequences with diversity
    - generate_genotype_matrix: Generate diploid/haploid genotype matrices
    - create_ecosystem: Create ecosystem with typed agents
    - run_simulation: Run ecosystem simulation steps
    - run_simulation_workflow: Execute configured simulation pipelines

Submodules:
    - models.sequences: DNA/protein sequence simulation
    - models.rna: RNA-seq count simulation
    - models.popgen: Population genetics simulations
    - models.agents: Agent-based modeling
    - workflow: Simulation workflow management
    - visualization: Plotting simulation results

Example:
    >>> from metainformant.simulation import generate_random_dna, mutate_sequence
    >>> seq = generate_random_dna(length=100)
    >>> mutated = mutate_sequence(seq, n_mut=5)

    >>> from metainformant.simulation import simulate_rnaseq_counts
    >>> counts = simulate_rnaseq_counts(n_genes=1000, n_samples=10)

    >>> from metainformant.simulation import create_ecosystem, run_simulation
    >>> eco = create_ecosystem(n_agents=50, agent_types=["producer", "consumer"], environment_size=(20, 20))
    >>> snapshots = run_simulation(eco, n_steps=100)

See Also:
    - docs/simulation/ for detailed simulation documentation
    - metainformant.math for population genetics theory
"""

from __future__ import annotations

# Import subpackages
from . import benchmark
from . import models
from . import visualization
from . import workflow

# Import modules from subpackages for backward compatibility
from .models import (
    agents,
    popgen,
    rna,
    sequences,
)

# Sequence simulation exports
from .models.sequences import (
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
    DNA_BASES,
    RNA_BASES,
    AMINO_ACIDS,
    GENETIC_CODE,
)

# RNA simulation exports
from .models.rna import (
    simulate_counts_negative_binomial,
    simulate_rnaseq_counts,
    simulate_differential_expression,
    simulate_bulk_rnaseq,
    simulate_single_cell_rnaseq,
    simulate_time_series_expression,
    simulate_spatial_expression,
    add_technical_noise,
)

# Population genetics exports
from .models.popgen import (
    generate_population_sequences,
    generate_two_populations,
    generate_genotype_matrix,
    generate_linkage_disequilibrium_data,
    simulate_bottleneck_population,
    simulate_population_expansion,
    generate_site_frequency_spectrum,
    simulate_admixture,
    simulate_selection,
)

# Agent simulation exports
from .models.agents import (
    Agent,
    GridAgent,
    GridWorld,
    Ecosystem,
    create_ecosystem,
    run_simulation,
    simulation_step,
    add_agent,
    remove_agent,
    get_population_dynamics,
    calculate_biodiversity_metrics,
    count_agents_by_type,
    simulate_predator_prey,
    simulate_competition,
)

# Workflow exports
from .workflow.workflow import (
    SimulationConfig,
    create_simulation_config,
    run_simulation_workflow,
    run_benchmark_simulation,
    validate_simulation_output,
    calibrate_simulation_parameters,
)

# Visualization exports
# Benchmark exports
from .benchmark.generators import (
    generate_benchmark_dataset,
    generate_synthetic_variants,
    generate_synthetic_expression,
    evaluate_benchmark,
    benchmark_suite,
)

from .visualization.visualization import (
    plot_sequence_evolution,
    animate_sequence_evolution,
    plot_rnaseq_simulation_results,
    plot_population_dynamics_simulation,
    plot_agent_based_model_results,
    plot_evolutionary_simulation_summary,
    plot_simulation_parameter_sensitivity,
    animate_population_dynamics,
    plot_simulation_validation_comparison,
    create_interactive_simulation_dashboard,
)

__all__ = [
    # Submodules
    "models",
    "visualization",
    "workflow",
    "agents",
    "popgen",
    "rna",
    "sequences",
    # Sequence simulation
    "generate_random_dna",
    "generate_random_protein",
    "mutate_sequence",
    "evolve_sequence",
    "translate_dna_to_protein",
    "reverse_transcribe_protein_to_dna",
    "generate_coding_sequence",
    "calculate_sequence_similarity",
    "generate_sequence_family",
    "analyze_sequence_divergence",
    "simulate_gene_duplication",
    "DNA_BASES",
    "RNA_BASES",
    "AMINO_ACIDS",
    "GENETIC_CODE",
    # RNA simulation
    "simulate_counts_negative_binomial",
    "simulate_rnaseq_counts",
    "simulate_differential_expression",
    "simulate_bulk_rnaseq",
    "simulate_single_cell_rnaseq",
    "simulate_time_series_expression",
    "simulate_spatial_expression",
    "add_technical_noise",
    # Population genetics
    "generate_population_sequences",
    "generate_two_populations",
    "generate_genotype_matrix",
    "generate_linkage_disequilibrium_data",
    "simulate_bottleneck_population",
    "simulate_population_expansion",
    "generate_site_frequency_spectrum",
    "simulate_admixture",
    "simulate_selection",
    # Agent simulation
    "Agent",
    "GridAgent",
    "GridWorld",
    "Ecosystem",
    "create_ecosystem",
    "run_simulation",
    "simulation_step",
    "add_agent",
    "remove_agent",
    "get_population_dynamics",
    "calculate_biodiversity_metrics",
    "count_agents_by_type",
    "simulate_predator_prey",
    "simulate_competition",
    # Workflow
    "SimulationConfig",
    "create_simulation_config",
    "run_simulation_workflow",
    "run_benchmark_simulation",
    "validate_simulation_output",
    "calibrate_simulation_parameters",
    # Visualization
    "plot_sequence_evolution",
    "animate_sequence_evolution",
    "plot_rnaseq_simulation_results",
    "plot_population_dynamics_simulation",
    "plot_agent_based_model_results",
    "plot_evolutionary_simulation_summary",
    "plot_simulation_parameter_sensitivity",
    "animate_population_dynamics",
    "plot_simulation_validation_comparison",
    "create_interactive_simulation_dashboard",
    # Benchmark
    "benchmark",
    "generate_benchmark_dataset",
    "generate_synthetic_variants",
    "generate_synthetic_expression",
    "evaluate_benchmark",
    "benchmark_suite",
]
