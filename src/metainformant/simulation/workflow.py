"""Simulation workflow orchestration and configuration management.

This module provides unified interfaces for running different types of
biological simulations with consistent configuration and validation.
All workflows support reproducible results and comprehensive output tracking.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Any, Callable
from dataclasses import dataclass, field
from datetime import datetime
import json

from metainformant.core import logging, errors, validation, config, io, paths

logger = logging.get_logger(__name__)


@dataclass
class SimulationConfig:
    """Configuration for simulation workflows."""

    # Simulation type
    simulation_type: str = "sequence_evolution"

    # Output settings
    output_dir: str | Path = "output/simulation"
    save_snapshots: bool = True
    snapshot_interval: int = 10

    # Random seed for reproducibility
    random_seed: Optional[int] = None

    # Common simulation parameters
    n_steps: int = 100
    population_size: int = 100
    mutation_rate: float = 0.001

    # Sequence simulation parameters
    sequence_length: int = 1000
    gc_content: float = 0.5

    # RNA simulation parameters
    n_genes: int = 1000
    n_samples: int = 50
    dispersion_mean: float = 0.5

    # Population genetics parameters
    n_snps: int = 1000
    selection_coefficient: float = 0.0

    # Agent-based simulation parameters
    n_agents: int = 100
    agent_types: List[str] = field(default_factory=lambda: ["producer", "consumer", "decomposer"])
    environment_size: Tuple[int, int] = (50, 50)

    # Validation and quality control
    validate_output: bool = True
    quality_checks: List[str] = field(default_factory=lambda: ["basic", "consistency"])

    def __post_init__(self):
        """Validate configuration after initialization."""
        valid_types = ["sequence_evolution", "population_genetics", "rna_expression",
                      "agent_ecosystem", "predator_prey", "competition"]
        if self.simulation_type not in valid_types:
            raise errors.ValidationError(f"Invalid simulation_type: {self.simulation_type}")

        validation.validate_range(self.n_steps, min_val=1, name="n_steps")
        validation.validate_range(self.population_size, min_val=2, name="population_size")
        validation.validate_range(self.mutation_rate, min_val=0.0, max_val=1.0, name="mutation_rate")
        validation.validate_range(self.sequence_length, min_val=1, name="sequence_length")
        validation.validate_range(self.gc_content, min_val=0.0, max_val=1.0, name="gc_content")
        validation.validate_range(self.n_genes, min_val=1, name="n_genes")
        validation.validate_range(self.n_samples, min_val=1, name="n_samples")
        validation.validate_range(self.dispersion_mean, min_val=0.0, name="dispersion_mean")
        validation.validate_range(self.n_snps, min_val=1, name="n_snps")
        validation.validate_range(self.n_agents, min_val=1, name="n_agents")
        validation.validate_type(self.agent_types, list, "agent_types")
        validation.validate_range(len(self.agent_types), min_val=1, name="agent_types length")


def create_simulation_config(simulation_type: str, parameters: Dict[str, Any]) -> SimulationConfig:
    """Create a simulation configuration with specified parameters.

    Args:
        simulation_type: Type of simulation to configure
        parameters: Dictionary of configuration parameters

    Returns:
        SimulationConfig object

    Raises:
        ValueError: If simulation_type is invalid
    """
    validation.validate_type(parameters, dict, "parameters")

    # Start with defaults for the simulation type
    config_dict = {
        "simulation_type": simulation_type,
        **parameters
    }

    try:
        config = SimulationConfig(**config_dict)
        logger.info(f"Created simulation config for type: {simulation_type}")
        return config
    except Exception as e:
        logger.error(f"Failed to create simulation config: {e}")
        raise errors.ConfigurationError(f"Invalid simulation configuration: {e}") from e


def run_benchmark_simulation(config: SimulationConfig, n_replicates: int = 10, *,
                           output_dir: str | Path | None = None) -> Dict[str, Any]:
    """Run multiple replicates of a simulation for benchmarking.

    Args:
        config: Simulation configuration
        n_replicates: Number of replicate runs
        output_dir: Directory for benchmark results

    Returns:
        Dictionary with benchmark results and statistics

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_replicates, min_val=1, name="n_replicates")

    if output_dir is None:
        output_dir = Path(config.output_dir) / "benchmark"
    else:
        output_dir = Path(output_dir)

    paths.ensure_directory(output_dir)

    logger.info(f"Starting benchmark simulation with {n_replicates} replicates")

    results = {
        "simulation_type": config.simulation_type,
        "n_replicates": n_replicates,
        "config": config.__dict__,
        "replicates": [],
        "summary": {},
        "timestamp": datetime.now().isoformat(),
    }

    all_execution_times = []
    all_outputs = []

    for rep in range(n_replicates):
        logger.info(f"Running replicate {rep + 1}/{n_replicates}")

        try:
            # Modify random seed for each replicate
            rep_config = SimulationConfig(**config.__dict__)
            if config.random_seed is not None:
                rep_config.random_seed = config.random_seed + rep

            # Run simulation
            start_time = datetime.now()
            rep_result = run_simulation_workflow(rep_config)
            end_time = datetime.now()

            execution_time = (end_time - start_time).total_seconds()

            replicate_data = {
                "replicate_id": rep,
                "execution_time": execution_time,
                "result": rep_result,
                "seed": rep_config.random_seed,
            }

            results["replicates"].append(replicate_data)
            all_execution_times.append(execution_time)
            all_outputs.append(rep_result)

        except Exception as e:
            error_msg = f"Replicate {rep} failed: {e}"
            logger.error(error_msg)
            results["replicates"].append({
                "replicate_id": rep,
                "error": error_msg,
                "execution_time": None,
            })

    # Calculate summary statistics
    successful_runs = [r for r in results["replicates"] if "error" not in r]

    if successful_runs:
        results["summary"] = {
            "total_replicates": n_replicates,
            "successful_replicates": len(successful_runs),
            "failed_replicates": n_replicates - len(successful_runs),
            "mean_execution_time": sum(all_execution_times) / len(all_execution_times),
            "min_execution_time": min(all_execution_times),
            "max_execution_time": max(all_execution_times),
            "execution_time_std": np.std(all_execution_times) if len(all_execution_times) > 1 else 0,
        }

        # Save individual replicate results
        for rep_data in results["replicates"]:
            if "error" not in rep_data:
                rep_file = output_dir / f"replicate_{rep_data['replicate_id']}.json"
                io.dump_json(rep_data, rep_file)

    # Save benchmark summary
    summary_file = output_dir / "benchmark_summary.json"
    io.dump_json(results, summary_file)

    logger.info(f"Benchmark completed. {len(successful_runs)}/{n_replicates} replicates successful")
    return results


def validate_simulation_output(simulation_data: Any, validation_criteria: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """Validate simulation output against specified criteria.

    Args:
        simulation_data: The simulation output to validate
        validation_criteria: Dictionary of validation criteria

    Returns:
        Tuple of (is_valid, list_of_issues)

    Raises:
        TypeError: If inputs are invalid
    """
    validation.validate_type(validation_criteria, dict, "validation_criteria")

    issues = []
    is_valid = True

    try:
        # Basic structure validation
        if not isinstance(simulation_data, dict):
            issues.append("Simulation output must be a dictionary")
            is_valid = False
            return is_valid, issues

        # Required fields validation
        required_fields = validation_criteria.get("required_fields", [])
        for field in required_fields:
            if field not in simulation_data:
                issues.append(f"Missing required field: {field}")
                is_valid = False

        # Data type validation
        type_checks = validation_criteria.get("type_checks", {})
        for field, expected_type in type_checks.items():
            if field in simulation_data:
                actual_value = simulation_data[field]
                if not isinstance(actual_value, expected_type):
                    issues.append(f"Field {field} has wrong type: expected {expected_type.__name__}, got {type(actual_value).__name__}")
                    is_valid = False

        # Range validation
        range_checks = validation_criteria.get("range_checks", {})
        for field, (min_val, max_val) in range_checks.items():
            if field in simulation_data:
                value = simulation_data[field]
                if isinstance(value, (int, float)):
                    if not (min_val <= value <= max_val):
                        issues.append(f"Field {field} out of range: {value} not in [{min_val}, {max_val}]")
                        is_valid = False

        # Custom validation functions
        custom_checks = validation_criteria.get("custom_checks", [])
        for check_func in custom_checks:
            if callable(check_func):
                try:
                    check_result = check_func(simulation_data)
                    if not check_result[0]:  # (is_valid, message)
                        issues.append(check_result[1])
                        is_valid = False
                except Exception as e:
                    issues.append(f"Custom validation failed: {e}")
                    is_valid = False

        # Consistency checks
        consistency_checks = validation_criteria.get("consistency_checks", [])
        for check in consistency_checks:
            if check == "positive_lengths" and "sequence_length" in simulation_data:
                if simulation_data["sequence_length"] <= 0:
                    issues.append("Sequence length must be positive")
                    is_valid = False

            elif check == "valid_probabilities":
                prob_fields = ["gc_content", "mutation_rate"]
                for field in prob_fields:
                    if field in simulation_data:
                        val = simulation_data[field]
                        if not (0.0 <= val <= 1.0):
                            issues.append(f"Probability field {field} must be in [0, 1]: got {val}")
                            is_valid = False

    except Exception as e:
        issues.append(f"Validation failed with exception: {e}")
        is_valid = False

    return is_valid, issues


def calibrate_simulation_parameters(target_data: Any, parameter_ranges: Dict[str, Tuple[float, float]], *,
                                  n_iterations: int = 100, fitness_function: Optional[Callable] = None) -> Dict[str, float]:
    """Calibrate simulation parameters to match target data.

    Args:
        target_data: Target data to match
        parameter_ranges: Dictionary mapping parameter names to (min, max) ranges
        n_iterations: Number of calibration iterations
        fitness_function: Function to evaluate parameter fitness

    Returns:
        Dictionary with calibrated parameters

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_iterations, min_val=1, name="n_iterations")
    validation.validate_type(parameter_ranges, dict, "parameter_ranges")

    if fitness_function is None:
        # Default fitness function - minimize difference
        def default_fitness(params, target):
            # Simple placeholder - would need domain-specific implementation
            return 1.0  # Perfect match

        fitness_function = default_fitness

    logger.info(f"Starting parameter calibration with {n_iterations} iterations")

    best_params = {}
    best_fitness = float('-inf')

    import random
    rng = random.Random()

    for iteration in range(n_iterations):
        # Sample random parameters
        params = {}
        for param_name, (min_val, max_val) in parameter_ranges.items():
            params[param_name] = rng.uniform(min_val, max_val)

        # Evaluate fitness
        try:
            fitness = fitness_function(params, target_data)

            if fitness > best_fitness:
                best_fitness = fitness
                best_params = params.copy()

        except Exception as e:
            logger.debug(f"Parameter evaluation failed for iteration {iteration}: {e}")

        if (iteration + 1) % 10 == 0:
            logger.info(f"Calibration iteration {iteration + 1}/{n_iterations}, best fitness: {best_fitness}")

    logger.info(f"Parameter calibration completed. Best fitness: {best_fitness}")
    return best_params


def run_simulation_workflow(config: SimulationConfig) -> Dict[str, Any]:
    """Run a simulation workflow based on configuration.

    Args:
        config: Simulation configuration

    Returns:
        Dictionary with simulation results

    Raises:
        ValueError: If simulation type is unsupported
    """
    validation.validate_type(config, SimulationConfig, "config")

    logger.info(f"Starting simulation workflow: {config.simulation_type}")

    # Set up random seed for reproducibility
    import random
    if config.random_seed is not None:
        random.seed(config.random_seed)
        import numpy as np
        np.random.seed(config.random_seed)

    rng = random.Random(config.random_seed)

    # Route to appropriate simulation function
    if config.simulation_type == "sequence_evolution":
        result = _run_sequence_simulation(config, rng)

    elif config.simulation_type == "population_genetics":
        result = _run_population_genetics_simulation(config, rng)

    elif config.simulation_type == "rna_expression":
        result = _run_rna_simulation(config, rng)

    elif config.simulation_type == "agent_ecosystem":
        result = _run_agent_simulation(config, rng)

    elif config.simulation_type == "predator_prey":
        result = _run_predator_prey_simulation(config, rng)

    elif config.simulation_type == "competition":
        result = _run_competition_simulation(config, rng)

    else:
        raise errors.ConfigurationError(f"Unsupported simulation type: {config.simulation_type}")

    # Add metadata
    result.update({
        "simulation_type": config.simulation_type,
        "config": config.__dict__,
        "random_seed": config.random_seed,
        "timestamp": datetime.now().isoformat(),
    })

    # Validate output if requested
    if config.validate_output:
        is_valid, issues = validate_simulation_output(result, {"required_fields": ["simulation_type"]})
        result["validation"] = {
            "is_valid": is_valid,
            "issues": issues,
        }

        if not is_valid:
            logger.warning(f"Simulation output validation failed: {issues}")

    # Save results if output directory specified
    if hasattr(config, 'output_dir') and config.output_dir:
        output_dir = Path(config.output_dir)
        paths.ensure_directory(output_dir)

        output_file = output_dir / f"{config.simulation_type}_result.json"
        io.dump_json(result, output_file)

        result["output_file"] = str(output_file)
        logger.info(f"Simulation results saved to {output_file}")

    logger.info(f"Simulation workflow completed: {config.simulation_type}")
    return result


def _run_sequence_simulation(config: SimulationConfig, rng: random.Random) -> Dict[str, Any]:
    """Run sequence evolution simulation."""
    from .sequences import generate_random_dna, evolve_sequence, analyze_sequence_divergence

    # Generate initial population
    population = []
    for _ in range(config.population_size):
        seq = generate_random_dna(config.sequence_length, gc_content=config.gc_content, rng=rng)
        population.append(seq)

    # Evolve population
    evolved_population = []
    for seq in population:
        evolved = evolve_sequence(seq, config.n_steps, mutation_rate=config.mutation_rate, rng=rng)
        evolved_population.append(evolved)

    # Analyze results
    divergence_analysis = analyze_sequence_divergence(evolved_population)

    return {
        "initial_population_size": config.population_size,
        "sequence_length": config.sequence_length,
        "generations": config.n_steps,
        "final_population": evolved_population[:5],  # Save first 5 sequences
        "divergence_analysis": divergence_analysis,
    }


def _run_population_genetics_simulation(config: SimulationConfig, rng: random.Random) -> Dict[str, Any]:
    """Run population genetics simulation."""
    from .popgen import generate_genotype_matrix, simulate_selection

    # Generate initial genotypes
    genotypes = generate_genotype_matrix(
        config.population_size, config.n_snps,
        maf_min=0.05, maf_max=0.5, rng=rng
    )

    # Apply selection if specified
    if config.selection_coefficient > 0:
        fitness_effects = np.random.normal(1.0, config.selection_coefficient, (config.n_snps, 3))
        selection_result = simulate_selection(
            genotypes, fitness_effects, config.n_steps,
            selection_strength=config.selection_coefficient, rng=rng
        )

        return {
            "population_size": config.population_size,
            "n_snps": config.n_snps,
            "generations": config.n_steps,
            "selection_coefficient": config.selection_coefficient,
            "initial_genotypes": genotypes.tolist(),
            "final_genotypes": selection_result["final_genotypes"].tolist(),
            "allele_frequencies": selection_result["allele_frequencies"],
        }
    else:
        return {
            "population_size": config.population_size,
            "n_snps": config.n_snps,
            "genotypes": genotypes.tolist(),
            "note": "Neutral evolution (no selection applied)",
        }


def _run_rna_simulation(config: SimulationConfig, rng: random.Random) -> Dict[str, Any]:
    """Run RNA expression simulation."""
    from .rna import simulate_bulk_rnaseq

    # Simulate expression data
    expression_matrix = simulate_bulk_rnaseq(
        config.n_samples, config.n_genes, rng=rng
    )

    return {
        "n_samples": config.n_samples,
        "n_genes": config.n_genes,
        "expression_matrix_shape": expression_matrix.shape,
        "total_reads": int(np.sum(expression_matrix)),
        "mean_expression_per_gene": np.mean(expression_matrix, axis=0).tolist()[:10],  # First 10 genes
        "expression_matrix_sample": expression_matrix[:5, :10].tolist(),  # Sample of data
    }


def _run_agent_simulation(config: SimulationConfig, rng: random.Random) -> Dict[str, Any]:
    """Run agent-based ecosystem simulation."""
    from .agents import create_ecosystem, run_simulation, get_population_dynamics, calculate_biodiversity_metrics

    # Create ecosystem
    ecosystem = create_ecosystem(
        config.n_agents, config.agent_types,
        config.environment_size, rng=rng
    )

    # Run simulation
    snapshots = run_simulation(
        ecosystem, config.n_steps,
        record_interval=config.snapshot_interval, rng=rng
    )

    # Analyze results
    population_dynamics = get_population_dynamics(snapshots)
    biodiversity = calculate_biodiversity_metrics(snapshots)

    return {
        "n_agents": config.n_agents,
        "agent_types": config.agent_types,
        "environment_size": config.environment_size,
        "generations": config.n_steps,
        "n_snapshots": len(snapshots),
        "population_dynamics": population_dynamics,
        "biodiversity_metrics": biodiversity,
        "final_snapshot": snapshots[-1] if snapshots else None,
    }


def _run_predator_prey_simulation(config: SimulationConfig, rng: random.Random) -> Dict[str, Any]:
    """Run predator-prey ecosystem simulation."""
    from .agents import create_ecosystem, simulate_predator_prey

    # Create ecosystem with predator-prey setup
    ecosystem = create_ecosystem(
        config.n_agents, ["producer", "consumer"],  # Prey and predators
        config.environment_size, rng=rng
    )

    # Run predator-prey simulation
    snapshots = simulate_predator_prey(
        ecosystem, config.n_steps,
        predator_efficiency=0.7, prey_reproduction=0.15, rng=rng
    )

    return {
        "simulation_type": "predator_prey",
        "n_agents": config.n_agents,
        "environment_size": config.environment_size,
        "generations": config.n_steps,
        "n_snapshots": len(snapshots),
        "predator_trajectory": [s["n_predators"] for s in snapshots],
        "prey_trajectory": [s["n_prey"] for s in snapshots],
        "final_predators": snapshots[-1]["n_predators"] if snapshots else 0,
        "final_prey": snapshots[-1]["n_prey"] if snapshots else 0,
    }


def _run_competition_simulation(config: SimulationConfig, rng: random.Random) -> Dict[str, Any]:
    """Run competition simulation."""
    from .agents import create_ecosystem, simulate_competition

    # Create ecosystem
    ecosystem = create_ecosystem(
        config.n_agents, config.agent_types,
        config.environment_size, rng=rng
    )

    # Run competition simulation
    snapshots = simulate_competition(
        ecosystem, config.n_steps,
        competition_strength=0.6, rng=rng
    )

    return {
        "simulation_type": "competition",
        "n_agents": config.n_agents,
        "agent_types": config.agent_types,
        "environment_size": config.environment_size,
        "generations": config.n_steps,
        "n_snapshots": len(snapshots),
        "occupied_positions_trajectory": [s["occupied_positions"] for s in snapshots],
        "competition_events_trajectory": [s["competition_events"] for s in snapshots],
        "final_agent_counts": snapshots[-1]["agent_counts"] if snapshots else {},
    }

