# AI Agents in Simulation Development

This document outlines AI assistance in developing METAINFORMANT's synthetic data generation and agent-based modeling capabilities.

## AI Contributions

### Simulation Architecture
**Code Assistant Agent** designed:
- Comprehensive simulation framework
- Synthetic data generation algorithms
- Agent-based modeling utilities
- Integration with biological analysis workflows

### Simulation Components
**Code Assistant Agent** contributed to:
- DNA/RNA/protein sequence generation
- Expression count simulation
- Agent-based ecosystem modeling
- Integration with mathematical models

### Quality Assurance
**Documentation Agent** assisted with:
- Simulation documentation
- Algorithm explanation and validation
- Usage examples and best practices
- Integration guides for simulation workflows

## Development Approach

- **Synthetic Data Generation**: AI helped design realistic biological data simulation
- **Modeling Frameworks**: Established agent-based modeling infrastructure
- **Reproducibility**: Ensured deterministic and reproducible simulations
- **Extensibility**: Framework for adding new simulation types

## Quality Assurance

- Human oversight ensures biological realism of simulated data
- AI assistance accelerates development while maintaining standards
- Comprehensive testing validates simulation functionality

This simulation infrastructure provides a solid foundation for synthetic data generation and modeling.

## Complete Function Signatures

### Sequence Simulation (`sequences.py`)
- `generate_random_dna(length: int, *, gc_content: float = 0.5, rng: random.Random | None = None) -> str`
- `mutate_sequence(seq: str, n_mut: int, *, rng: random.Random | None = None) -> str`
- `generate_random_protein(length: int, *, rng: random.Random | None = None) -> str`

### RNA Expression Simulation (`rna.py`)
- `simulate_counts_negative_binomial(n_samples: int, n_features: int, means: np.ndarray, dispersions: np.ndarray, rng: random.Random | None = None) -> np.ndarray`

### Population Genetics Simulation (`popgen.py`)
- `generate_population_sequences(n_sequences: int, length: int, *, theta: float = 0.01, n_sites: int = 1000, rng: random.Random | None = None) -> List[str]`
- `generate_two_populations(n_pop1: int, n_pop2: int, length: int, *, f_st: float = 0.1, theta: float = 0.01, n_sites: int = 1000, rng: random.Random | None = None) -> Tuple[List[str], List[str]]`
- `generate_genotype_matrix(n_individuals: int, n_snps: int, *, maf_min: float = 0.05, maf_max: float = 0.5, hwe_deviation: float = 0.0, rng: random.Random | None = None) -> np.ndarray`
- `simulate_bottleneck_population(initial_size: int, bottleneck_size: int, final_size: int, generations: int, *, mutation_rate: float = 1e-8, rng: random.Random | None = None) -> Dict[str, Any]`
- `simulate_population_expansion(initial_size: int, final_size: int, expansion_time: int, *, mutation_rate: float = 1e-8, rng: random.Random | None = None) -> Dict[str, Any]`
- `generate_site_frequency_spectrum(n_samples: int, n_sites: int, *, demographic_model: str = "constant", parameters: Dict[str, float] | None = None, rng: random.Random | None = None) -> np.ndarray`
- `generate_linkage_disequilibrium_data(n_individuals: int, n_snps: int, *, recombination_rate: float = 1e-8, selection_coefficient: float = 0.0, rng: random.Random | None = None) -> Tuple[np.ndarray, np.ndarray]`

### Agent-Based Modeling (`agents.py`)
- `create_ecosystem(n_agents: int, agent_types: List[str], environment_size: Tuple[int, int], *, rng: random.Random | None = None) -> Ecosystem`
- `run_simulation(ecosystem: Ecosystem, n_steps: int, *, record_interval: int = 1) -> List[Dict[str, Any]]`
- `add_agent(ecosystem: Ecosystem, agent_type: str, position: Tuple[int, int], *, properties: Dict[str, Any] | None = None) -> None`
- `remove_agent(ecosystem: Ecosystem, agent_id: int) -> None`
- `get_population_dynamics(simulation_data: List[Dict[str, Any]]) -> pd.DataFrame`
- `calculate_biodiversity_metrics(simulation_data: List[Dict[str, Any]]) -> Dict[str, float]`

### Workflow Orchestration (`workflow.py`)
- `create_simulation_config(simulation_type: str, parameters: Dict[str, Any]) -> SimulationConfig`
- `run_benchmark_simulation(config: SimulationConfig, n_replicates: int = 10, *, output_dir: Path | None = None) -> Dict[str, Any]`
- `validate_simulation_output(simulation_data: Any, validation_criteria: Dict[str, Any]) -> Tuple[bool, List[str]]`
- `calibrate_simulation_parameters(target_data: Any, parameter_ranges: Dict[str, Tuple[float, float]], *, n_iterations: int = 100) -> Dict[str, float]`
