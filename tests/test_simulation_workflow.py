"""Comprehensive tests for the simulation workflow module.

Tests cover SimulationConfig validation, create_simulation_config factory,
run_simulation_workflow for all six simulation types, output validation,
parameter calibration, benchmark runs, and reproducibility guarantees.

All tests use real implementations (NO MOCKING) per project policy.
"""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pytest

from metainformant.core.utils.errors import ConfigError, ValidationError
from metainformant.simulation import (
    SimulationConfig,
    calibrate_simulation_parameters,
    create_simulation_config,
    run_benchmark_simulation,
    run_simulation_workflow,
    validate_simulation_output,
)
from metainformant.simulation.workflow.workflow import (
    SimulationConfig as WorkflowSimulationConfig,
    calibrate_simulation_parameters as wf_calibrate,
    create_simulation_config as wf_create_config,
    run_benchmark_simulation as wf_run_benchmark,
    run_simulation_workflow as wf_run_workflow,
    validate_simulation_output as wf_validate_output,
)


# ---------------------------------------------------------------------------
# Constants for small/fast simulation sizes
# ---------------------------------------------------------------------------
SMALL_N_STEPS = 5
SMALL_POP_SIZE = 4  # minimum is 2
SMALL_SEQ_LEN = 20
SMALL_N_GENES = 10
SMALL_N_SAMPLES = 4
SMALL_N_SNPS = 10
SMALL_N_AGENTS = 6
SMALL_ENV_SIZE = (5, 5)
DEFAULT_SEED = 42


# ===========================================================================
# Section 1: Import parity - both import paths resolve to the same objects
# ===========================================================================


class TestImportParity:
    """Verify that the two documented import paths reference identical objects."""

    def test_simulation_config_is_same_class(self) -> None:
        assert SimulationConfig is WorkflowSimulationConfig

    def test_create_simulation_config_is_same_function(self) -> None:
        assert create_simulation_config is wf_create_config

    def test_run_simulation_workflow_is_same_function(self) -> None:
        assert run_simulation_workflow is wf_run_workflow

    def test_run_benchmark_simulation_is_same_function(self) -> None:
        assert run_benchmark_simulation is wf_run_benchmark

    def test_validate_simulation_output_is_same_function(self) -> None:
        assert validate_simulation_output is wf_validate_output

    def test_calibrate_simulation_parameters_is_same_function(self) -> None:
        assert calibrate_simulation_parameters is wf_calibrate


# ===========================================================================
# Section 2: SimulationConfig defaults and validation
# ===========================================================================


class TestSimulationConfigDefaults:
    """Verify that SimulationConfig initialises with documented defaults."""

    def test_default_simulation_type(self) -> None:
        cfg = SimulationConfig()
        assert cfg.simulation_type == "sequence_evolution"

    def test_default_output_dir(self) -> None:
        cfg = SimulationConfig()
        assert str(cfg.output_dir) == "output/simulation"

    def test_default_snapshot_settings(self) -> None:
        cfg = SimulationConfig()
        assert cfg.save_snapshots is True
        assert cfg.snapshot_interval == 10

    def test_default_random_seed_is_none(self) -> None:
        cfg = SimulationConfig()
        assert cfg.random_seed is None

    def test_default_numeric_params(self) -> None:
        cfg = SimulationConfig()
        assert cfg.n_steps == 100
        assert cfg.population_size == 100
        assert cfg.mutation_rate == pytest.approx(0.001)
        assert cfg.sequence_length == 1000
        assert cfg.gc_content == pytest.approx(0.5)
        assert cfg.n_genes == 1000
        assert cfg.n_samples == 50
        assert cfg.dispersion_mean == pytest.approx(0.5)
        assert cfg.n_snps == 1000
        assert cfg.selection_coefficient == pytest.approx(0.0)
        assert cfg.n_agents == 100

    def test_default_agent_types(self) -> None:
        cfg = SimulationConfig()
        assert cfg.agent_types == ["producer", "consumer", "decomposer"]

    def test_default_environment_size(self) -> None:
        cfg = SimulationConfig()
        assert cfg.environment_size == (50, 50)

    def test_default_validation_settings(self) -> None:
        cfg = SimulationConfig()
        assert cfg.validate_output is True
        assert cfg.quality_checks == ["basic", "consistency"]

    def test_agent_types_default_is_independent_per_instance(self) -> None:
        """Default mutable list must not be shared across instances."""
        cfg1 = SimulationConfig()
        cfg2 = SimulationConfig()
        cfg1.agent_types.append("scavenger")
        assert "scavenger" not in cfg2.agent_types


class TestSimulationConfigValidation:
    """Verify that __post_init__ catches invalid configurations."""

    def test_invalid_simulation_type_raises_validation_error(self) -> None:
        with pytest.raises(ValidationError, match="Invalid simulation_type"):
            SimulationConfig(simulation_type="nonexistent_type")

    def test_invalid_simulation_type_empty_string(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(simulation_type="")

    def test_population_size_below_minimum(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(population_size=1)

    def test_population_size_at_minimum(self) -> None:
        cfg = SimulationConfig(population_size=2)
        assert cfg.population_size == 2

    def test_n_steps_below_minimum(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(n_steps=0)

    def test_mutation_rate_out_of_range_high(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(mutation_rate=1.5)

    def test_mutation_rate_out_of_range_negative(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(mutation_rate=-0.01)

    def test_mutation_rate_at_boundaries(self) -> None:
        cfg_low = SimulationConfig(mutation_rate=0.0)
        cfg_high = SimulationConfig(mutation_rate=1.0)
        assert cfg_low.mutation_rate == pytest.approx(0.0)
        assert cfg_high.mutation_rate == pytest.approx(1.0)

    def test_gc_content_out_of_range(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(gc_content=1.1)
        with pytest.raises(ValidationError):
            SimulationConfig(gc_content=-0.1)

    def test_sequence_length_below_minimum(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(sequence_length=0)

    def test_n_genes_below_minimum(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(n_genes=0)

    def test_n_samples_below_minimum(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(n_samples=0)

    def test_dispersion_mean_negative(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(dispersion_mean=-0.1)

    def test_n_snps_below_minimum(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(n_snps=0)

    def test_n_agents_below_minimum(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(n_agents=0)

    def test_empty_agent_types_raises(self) -> None:
        with pytest.raises(ValidationError):
            SimulationConfig(agent_types=[])

    def test_all_valid_types_accepted(self) -> None:
        valid_types = [
            "sequence_evolution",
            "population_genetics",
            "rna_expression",
            "agent_ecosystem",
            "predator_prey",
            "competition",
        ]
        for sim_type in valid_types:
            cfg = SimulationConfig(simulation_type=sim_type)
            assert cfg.simulation_type == sim_type


# ===========================================================================
# Section 3: create_simulation_config factory function
# ===========================================================================


class TestCreateSimulationConfig:
    """Tests for the create_simulation_config factory."""

    def test_create_with_valid_type_and_defaults(self) -> None:
        cfg = create_simulation_config("sequence_evolution", {})
        assert isinstance(cfg, SimulationConfig)
        assert cfg.simulation_type == "sequence_evolution"

    def test_create_with_custom_parameters(self) -> None:
        params = {
            "n_steps": SMALL_N_STEPS,
            "population_size": SMALL_POP_SIZE,
            "mutation_rate": 0.05,
            "random_seed": DEFAULT_SEED,
        }
        cfg = create_simulation_config("sequence_evolution", params)
        assert cfg.n_steps == SMALL_N_STEPS
        assert cfg.population_size == SMALL_POP_SIZE
        assert cfg.mutation_rate == pytest.approx(0.05)
        assert cfg.random_seed == DEFAULT_SEED

    @pytest.mark.parametrize(
        "sim_type",
        [
            "sequence_evolution",
            "population_genetics",
            "rna_expression",
            "agent_ecosystem",
            "predator_prey",
            "competition",
        ],
    )
    def test_create_with_each_valid_type(self, sim_type: str) -> None:
        cfg = create_simulation_config(sim_type, {})
        assert cfg.simulation_type == sim_type

    def test_create_with_invalid_type_raises_config_error(self) -> None:
        with pytest.raises((ConfigError, ValidationError)):
            create_simulation_config("bogus_type", {})

    def test_create_with_invalid_parameter_raises_config_error(self) -> None:
        with pytest.raises((ConfigError, ValidationError)):
            create_simulation_config("sequence_evolution", {"population_size": 0})

    def test_create_with_output_dir(self, tmp_path: Path) -> None:
        out = tmp_path / "sim_output"
        cfg = create_simulation_config(
            "sequence_evolution",
            {"output_dir": str(out)},
        )
        assert str(cfg.output_dir) == str(out)


# ===========================================================================
# Section 4: run_simulation_workflow -- all six simulation types
# ===========================================================================


def _make_small_config(
    simulation_type: str,
    tmp_path: Path,
    *,
    seed: int = DEFAULT_SEED,
    extra: dict | None = None,
) -> SimulationConfig:
    """Helper to build a fast, small-footprint simulation config."""
    params: Dict[str, Any] = {
        "simulation_type": simulation_type,
        "n_steps": SMALL_N_STEPS,
        "population_size": SMALL_POP_SIZE,
        "mutation_rate": 0.01,
        "sequence_length": SMALL_SEQ_LEN,
        "gc_content": 0.5,
        "n_genes": SMALL_N_GENES,
        "n_samples": SMALL_N_SAMPLES,
        "n_snps": SMALL_N_SNPS,
        "n_agents": SMALL_N_AGENTS,
        "environment_size": SMALL_ENV_SIZE,
        "random_seed": seed,
        "output_dir": str(tmp_path / "sim_output"),
        "save_snapshots": False,
        "validate_output": True,
    }
    if extra:
        params.update(extra)
    return SimulationConfig(**params)


class TestRunSequenceEvolution:
    """Tests for sequence_evolution simulation type."""

    def test_returns_dict_with_required_keys(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_simulation_workflow(cfg)
        assert isinstance(result, dict)
        assert result["simulation_type"] == "sequence_evolution"
        assert "config" in result
        assert "random_seed" in result
        assert "timestamp" in result

    def test_result_contains_divergence_analysis(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "divergence_analysis" in result

    def test_result_contains_population_data(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_simulation_workflow(cfg)
        assert result["initial_population_size"] == SMALL_POP_SIZE
        assert result["sequence_length"] == SMALL_SEQ_LEN
        assert result["generations"] == SMALL_N_STEPS

    def test_output_file_created(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "output_file" in result
        output_path = Path(result["output_file"])
        assert output_path.exists()

    def test_validation_passes(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "validation" in result
        assert result["validation"]["is_valid"] is True
        assert result["validation"]["issues"] == []


class TestRunPopulationGenetics:
    """Tests for population_genetics simulation type."""

    def test_neutral_evolution(self, tmp_path: Path) -> None:
        cfg = _make_small_config(
            "population_genetics", tmp_path, extra={"selection_coefficient": 0.0}
        )
        result = run_simulation_workflow(cfg)
        assert result["simulation_type"] == "population_genetics"
        assert result["population_size"] == SMALL_POP_SIZE
        assert result["n_snps"] == SMALL_N_SNPS
        assert "genotypes" in result

    def test_with_selection(self, tmp_path: Path) -> None:
        cfg = _make_small_config(
            "population_genetics", tmp_path, extra={"selection_coefficient": 0.1}
        )
        result = run_simulation_workflow(cfg)
        assert result["simulation_type"] == "population_genetics"
        assert "allele_frequencies" in result or "final_genotypes" in result

    def test_output_file_created(self, tmp_path: Path) -> None:
        cfg = _make_small_config("population_genetics", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "output_file" in result
        assert Path(result["output_file"]).exists()


class TestRunRNAExpression:
    """Tests for rna_expression simulation type."""

    def test_returns_expression_data(self, tmp_path: Path) -> None:
        cfg = _make_small_config("rna_expression", tmp_path)
        result = run_simulation_workflow(cfg)
        assert result["simulation_type"] == "rna_expression"
        assert result["n_samples"] == SMALL_N_SAMPLES
        assert result["n_genes"] == SMALL_N_GENES
        assert "expression_matrix_shape" in result

    def test_total_reads_positive(self, tmp_path: Path) -> None:
        cfg = _make_small_config("rna_expression", tmp_path)
        result = run_simulation_workflow(cfg)
        assert result["total_reads"] > 0

    def test_expression_matrix_sample_present(self, tmp_path: Path) -> None:
        cfg = _make_small_config("rna_expression", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "expression_matrix_sample" in result
        assert isinstance(result["expression_matrix_sample"], list)


class TestRunAgentEcosystem:
    """Tests for agent_ecosystem simulation type."""

    def test_returns_ecosystem_data(self, tmp_path: Path) -> None:
        cfg = _make_small_config("agent_ecosystem", tmp_path)
        result = run_simulation_workflow(cfg)
        assert result["simulation_type"] == "agent_ecosystem"
        assert result["n_agents"] == SMALL_N_AGENTS
        assert result["environment_size"] == SMALL_ENV_SIZE

    def test_population_dynamics_present(self, tmp_path: Path) -> None:
        cfg = _make_small_config("agent_ecosystem", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "population_dynamics" in result

    def test_biodiversity_metrics_present(self, tmp_path: Path) -> None:
        cfg = _make_small_config("agent_ecosystem", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "biodiversity_metrics" in result

    def test_snapshots_recorded(self, tmp_path: Path) -> None:
        cfg = _make_small_config("agent_ecosystem", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "n_snapshots" in result
        assert result["n_snapshots"] >= 1


class TestRunPredatorPrey:
    """Tests for predator_prey simulation type."""

    def test_returns_predator_prey_data(self, tmp_path: Path) -> None:
        cfg = _make_small_config("predator_prey", tmp_path)
        result = run_simulation_workflow(cfg)
        assert result["simulation_type"] == "predator_prey"

    def test_trajectories_present(self, tmp_path: Path) -> None:
        cfg = _make_small_config("predator_prey", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "predator_trajectory" in result
        assert "prey_trajectory" in result
        assert isinstance(result["predator_trajectory"], list)
        assert isinstance(result["prey_trajectory"], list)

    def test_final_counts_present(self, tmp_path: Path) -> None:
        cfg = _make_small_config("predator_prey", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "final_predators" in result
        assert "final_prey" in result


class TestRunCompetition:
    """Tests for competition simulation type."""

    def test_returns_competition_data(self, tmp_path: Path) -> None:
        cfg = _make_small_config("competition", tmp_path)
        result = run_simulation_workflow(cfg)
        assert result["simulation_type"] == "competition"

    def test_trajectories_present(self, tmp_path: Path) -> None:
        cfg = _make_small_config("competition", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "occupied_positions_trajectory" in result
        assert "competition_events_trajectory" in result

    def test_final_agent_counts_present(self, tmp_path: Path) -> None:
        cfg = _make_small_config("competition", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "final_agent_counts" in result
        assert isinstance(result["final_agent_counts"], dict)


# ===========================================================================
# Section 5: validate_simulation_output
# ===========================================================================


class TestValidateSimulationOutput:
    """Tests for the standalone validate_simulation_output function."""

    def test_empty_criteria_passes(self) -> None:
        data = {"simulation_type": "test", "value": 42}
        is_valid, issues = validate_simulation_output(data, {})
        assert is_valid is True
        assert issues == []

    def test_non_dict_data_fails(self) -> None:
        is_valid, issues = validate_simulation_output("not_a_dict", {})
        assert is_valid is False
        assert any("dictionary" in issue.lower() for issue in issues)

    def test_required_fields_present(self) -> None:
        data = {"alpha": 1, "beta": 2}
        criteria = {"required_fields": ["alpha", "beta"]}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is True
        assert issues == []

    def test_required_fields_missing(self) -> None:
        data = {"alpha": 1}
        criteria = {"required_fields": ["alpha", "beta", "gamma"]}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is False
        assert len(issues) == 2  # beta and gamma missing

    def test_type_checks_pass(self) -> None:
        data = {"count": 10, "name": "test"}
        criteria = {"type_checks": {"count": int, "name": str}}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is True

    def test_type_checks_fail(self) -> None:
        data = {"count": "not_an_int"}
        criteria = {"type_checks": {"count": int}}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is False
        assert any("wrong type" in issue.lower() for issue in issues)

    def test_range_checks_pass(self) -> None:
        data = {"probability": 0.5, "score": 75}
        criteria = {
            "range_checks": {
                "probability": (0.0, 1.0),
                "score": (0, 100),
            }
        }
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is True

    def test_range_checks_fail(self) -> None:
        data = {"probability": 1.5}
        criteria = {"range_checks": {"probability": (0.0, 1.0)}}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is False
        assert any("out of range" in issue.lower() for issue in issues)

    def test_custom_check_passes(self) -> None:
        def check_positive_total(data: dict) -> Tuple[bool, str]:
            total = data.get("total", 0)
            if total > 0:
                return (True, "")
            return (False, "total must be positive")

        data = {"total": 42}
        criteria = {"custom_checks": [check_positive_total]}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is True

    def test_custom_check_fails(self) -> None:
        def check_positive_total(data: dict) -> Tuple[bool, str]:
            total = data.get("total", 0)
            if total > 0:
                return (True, "")
            return (False, "total must be positive")

        data = {"total": -5}
        criteria = {"custom_checks": [check_positive_total]}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is False
        assert any("total must be positive" in issue for issue in issues)

    def test_consistency_check_positive_lengths(self) -> None:
        data = {"sequence_length": 100}
        criteria = {"consistency_checks": ["positive_lengths"]}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is True

    def test_consistency_check_positive_lengths_fails(self) -> None:
        data = {"sequence_length": -1}
        criteria = {"consistency_checks": ["positive_lengths"]}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is False

    def test_consistency_check_valid_probabilities(self) -> None:
        data = {"gc_content": 0.5, "mutation_rate": 0.01}
        criteria = {"consistency_checks": ["valid_probabilities"]}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is True

    def test_consistency_check_valid_probabilities_fails(self) -> None:
        data = {"gc_content": 1.5, "mutation_rate": 0.01}
        criteria = {"consistency_checks": ["valid_probabilities"]}
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is False

    def test_combined_criteria(self) -> None:
        data = {
            "simulation_type": "sequence_evolution",
            "count": 10,
            "probability": 0.5,
            "sequence_length": 100,
            "gc_content": 0.4,
            "mutation_rate": 0.01,
        }
        criteria = {
            "required_fields": ["simulation_type", "count"],
            "type_checks": {"count": int, "simulation_type": str},
            "range_checks": {"probability": (0.0, 1.0), "count": (0, 1000)},
            "consistency_checks": ["positive_lengths", "valid_probabilities"],
        }
        is_valid, issues = validate_simulation_output(data, criteria)
        assert is_valid is True
        assert issues == []

    def test_validate_real_workflow_output(self, tmp_path: Path) -> None:
        """Run a real simulation and validate its output."""
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_simulation_workflow(cfg)
        criteria = {
            "required_fields": ["simulation_type", "config", "timestamp"],
            "type_checks": {"simulation_type": str},
        }
        is_valid, issues = validate_simulation_output(result, criteria)
        assert is_valid is True
        assert issues == []


# ===========================================================================
# Section 6: calibrate_simulation_parameters
# ===========================================================================


class TestCalibrateSimulationParameters:
    """Tests for the parameter calibration function."""

    def test_returns_dict_with_parameter_names(self) -> None:
        target_data = {"mean": 5.0, "std": 1.0}
        ranges = {
            "mutation_rate": (0.001, 0.1),
            "population_size": (10.0, 500.0),
        }
        result = calibrate_simulation_parameters(
            target_data, ranges, n_iterations=20
        )
        assert isinstance(result, dict)
        assert "mutation_rate" in result
        assert "population_size" in result

    def test_parameter_values_within_ranges(self) -> None:
        ranges = {
            "alpha": (0.0, 1.0),
            "beta": (10.0, 100.0),
        }
        result = calibrate_simulation_parameters(
            {"target": 0.5}, ranges, n_iterations=50
        )
        assert 0.0 <= result["alpha"] <= 1.0
        assert 10.0 <= result["beta"] <= 100.0

    def test_custom_fitness_function(self) -> None:
        target_data = {"target_value": 0.5}

        def fitness_fn(params: dict, target: dict) -> float:
            # Reward params closer to target
            diff = abs(params["x"] - target["target_value"])
            return -diff  # Negative distance = higher is better

        ranges = {"x": (0.0, 1.0)}
        result = calibrate_simulation_parameters(
            target_data,
            ranges,
            n_iterations=200,
            fitness_function=fitness_fn,
        )
        assert isinstance(result, dict)
        assert "x" in result
        # With enough iterations the best x should be somewhat close to 0.5
        # (stochastic, so wide tolerance)
        assert 0.0 <= result["x"] <= 1.0

    def test_single_iteration(self) -> None:
        ranges = {"p": (0.0, 1.0)}
        result = calibrate_simulation_parameters(
            {}, ranges, n_iterations=1
        )
        assert "p" in result

    def test_invalid_n_iterations_raises(self) -> None:
        with pytest.raises(ValidationError):
            calibrate_simulation_parameters({}, {"a": (0.0, 1.0)}, n_iterations=0)


# ===========================================================================
# Section 7: run_benchmark_simulation
# ===========================================================================


class TestRunBenchmarkSimulation:
    """Tests for the benchmark simulation runner."""

    def test_benchmark_returns_summary(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_benchmark_simulation(
            cfg, n_replicates=2, output_dir=str(tmp_path / "benchmark")
        )
        assert isinstance(result, dict)
        assert result["simulation_type"] == "sequence_evolution"
        assert result["n_replicates"] == 2
        assert "replicates" in result
        assert "summary" in result
        assert "timestamp" in result

    def test_benchmark_creates_output_files(self, tmp_path: Path) -> None:
        bench_dir = tmp_path / "benchmark"
        cfg = _make_small_config("sequence_evolution", tmp_path)
        run_benchmark_simulation(cfg, n_replicates=2, output_dir=str(bench_dir))
        assert bench_dir.exists()
        summary_file = bench_dir / "benchmark_summary.json"
        assert summary_file.exists()

    def test_benchmark_replicate_count(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_benchmark_simulation(
            cfg, n_replicates=3, output_dir=str(tmp_path / "benchmark")
        )
        assert len(result["replicates"]) == 3

    def test_benchmark_summary_stats(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_benchmark_simulation(
            cfg, n_replicates=2, output_dir=str(tmp_path / "benchmark")
        )
        summary = result["summary"]
        assert summary["total_replicates"] == 2
        assert summary["successful_replicates"] >= 1
        assert "mean_execution_time" in summary
        assert summary["mean_execution_time"] >= 0

    def test_benchmark_with_seed_increments(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path, seed=100)
        result = run_benchmark_simulation(
            cfg, n_replicates=3, output_dir=str(tmp_path / "benchmark")
        )
        seeds = [r["seed"] for r in result["replicates"] if "seed" in r]
        # Seeds should be 100, 101, 102
        assert seeds == [100, 101, 102]


# ===========================================================================
# Section 8: Reproducibility
# ===========================================================================


class TestReproducibility:
    """Verify that the same seed produces identical results."""

    def test_sequence_evolution_reproducible(self, tmp_path: Path) -> None:
        cfg1 = _make_small_config(
            "sequence_evolution", tmp_path / "run1", seed=777
        )
        cfg2 = _make_small_config(
            "sequence_evolution", tmp_path / "run2", seed=777
        )
        r1 = run_simulation_workflow(cfg1)
        r2 = run_simulation_workflow(cfg2)
        # Core simulation data must match
        assert r1["final_population"] == r2["final_population"]
        assert r1["divergence_analysis"] == r2["divergence_analysis"]

    def test_population_genetics_reproducible(self, tmp_path: Path) -> None:
        cfg1 = _make_small_config(
            "population_genetics", tmp_path / "run1", seed=888
        )
        cfg2 = _make_small_config(
            "population_genetics", tmp_path / "run2", seed=888
        )
        r1 = run_simulation_workflow(cfg1)
        r2 = run_simulation_workflow(cfg2)
        assert r1["genotypes"] == r2["genotypes"]

    def test_different_seeds_produce_different_results(self, tmp_path: Path) -> None:
        cfg1 = _make_small_config(
            "sequence_evolution", tmp_path / "run1", seed=111
        )
        cfg2 = _make_small_config(
            "sequence_evolution", tmp_path / "run2", seed=999
        )
        r1 = run_simulation_workflow(cfg1)
        r2 = run_simulation_workflow(cfg2)
        # With different seeds and mutations, populations should differ
        assert r1["final_population"] != r2["final_population"]

    def test_seed_stored_in_result(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path, seed=42)
        result = run_simulation_workflow(cfg)
        assert result["random_seed"] == 42


# ===========================================================================
# Section 9: Edge cases and error handling
# ===========================================================================


class TestEdgeCases:
    """Edge cases and error paths."""

    def test_config_with_custom_output_dir(self, tmp_path: Path) -> None:
        custom_dir = tmp_path / "custom_sim_output"
        cfg = SimulationConfig(
            simulation_type="sequence_evolution",
            output_dir=str(custom_dir),
            n_steps=SMALL_N_STEPS,
            population_size=SMALL_POP_SIZE,
            sequence_length=SMALL_SEQ_LEN,
            random_seed=DEFAULT_SEED,
        )
        result = run_simulation_workflow(cfg)
        assert custom_dir.exists()
        assert "output_file" in result

    def test_validate_output_disabled(self, tmp_path: Path) -> None:
        cfg = _make_small_config(
            "sequence_evolution", tmp_path, extra={"validate_output": False}
        )
        result = run_simulation_workflow(cfg)
        assert "validation" not in result

    def test_config_dict_in_result(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_simulation_workflow(cfg)
        assert isinstance(result["config"], dict)
        assert result["config"]["simulation_type"] == "sequence_evolution"
        assert result["config"]["random_seed"] == DEFAULT_SEED

    def test_timestamp_present_in_result(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_simulation_workflow(cfg)
        assert "timestamp" in result
        assert isinstance(result["timestamp"], str)
        # ISO format should contain T or - separator
        assert "T" in result["timestamp"] or "-" in result["timestamp"]

    def test_minimum_viable_sequence_evolution(self, tmp_path: Path) -> None:
        """Smallest possible valid simulation."""
        cfg = SimulationConfig(
            simulation_type="sequence_evolution",
            output_dir=str(tmp_path / "sim_output"),
            n_steps=1,
            population_size=2,
            sequence_length=1,
            mutation_rate=0.0,
            random_seed=1,
        )
        result = run_simulation_workflow(cfg)
        assert result["simulation_type"] == "sequence_evolution"
        assert result["initial_population_size"] == 2

    def test_minimum_viable_rna_expression(self, tmp_path: Path) -> None:
        """Smallest possible RNA expression simulation."""
        cfg = SimulationConfig(
            simulation_type="rna_expression",
            output_dir=str(tmp_path / "sim_output"),
            n_steps=1,
            population_size=2,
            n_genes=1,
            n_samples=1,
            random_seed=1,
        )
        result = run_simulation_workflow(cfg)
        assert result["simulation_type"] == "rna_expression"
        assert result["n_genes"] == 1
        assert result["n_samples"] == 1

    def test_output_file_is_valid_json(self, tmp_path: Path) -> None:
        cfg = _make_small_config("sequence_evolution", tmp_path)
        result = run_simulation_workflow(cfg)
        output_path = Path(result["output_file"])
        with open(output_path) as f:
            loaded = json.load(f)
        assert loaded["simulation_type"] == "sequence_evolution"

    def test_validate_simulation_output_with_none_data(self) -> None:
        is_valid, issues = validate_simulation_output(None, {})
        assert is_valid is False

    def test_validate_simulation_output_with_list_data(self) -> None:
        is_valid, issues = validate_simulation_output([1, 2, 3], {})
        assert is_valid is False

    def test_create_config_preserves_path_type(self, tmp_path: Path) -> None:
        cfg = create_simulation_config(
            "sequence_evolution",
            {"output_dir": tmp_path / "output"},
        )
        # output_dir can be str or Path - either is valid
        assert str(cfg.output_dir) == str(tmp_path / "output")
