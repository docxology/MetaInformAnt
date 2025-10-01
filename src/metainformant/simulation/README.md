# Simulation Module

The `simulation` module provides synthetic data generation and agent-based modeling tools for testing hypotheses and generating training data across biological domains.

## Overview

This module enables researchers to:
- **Generate Synthetic Data**: Create realistic biological sequences and expression data
- **Agent-Based Modeling**: Simulate complex biological systems with interacting agents
- **Hypothesis Testing**: Test theoretical models with controlled synthetic data
- **Training Data Generation**: Create datasets for machine learning applications

## Submodules

### Sequence Simulation (`sequences.py`)
Synthetic DNA, RNA, and protein sequence generation.

**Key Features:**
- Random sequence generation with biological constraints
- Mutation simulation with configurable rates
- Sequence evolution under selection
- Realistic length distributions and composition

**Usage:**
```python
from metainformant.simulation import sequences

# Generate random sequences
dna_seq = sequences.generate_random_dna(1000)
rna_seq = sequences.generate_random_rna(500)
protein_seq = sequences.generate_random_protein(200)

# Simulate mutations
mutated = sequences.mutate_sequence(dna_seq, mutation_rate=0.01)
```

### Expression Simulation (`rna.py`)
Synthetic gene expression data generation.

**Key Features:**
- Negative binomial distribution modeling
- Library size and dispersion control
- Batch effect simulation
- Differential expression simulation

**Usage:**
```python
from metainformant.simulation import rna

# Simulate expression counts
counts = rna.simulate_counts_negative_binomial(
    n_genes=1000,
    n_samples=20,
    library_size=1000000
)

# Add differential expression
diff_counts = rna.simulate_differential_expression(
    counts,
    n_de_genes=50,
    fold_changes=[0.5, 2.0]
)
```

### Agent-Based Models (`agents.py`)
Grid-world and multi-agent simulation framework.

**Key Features:**
- Configurable grid environments
- Agent behavior modeling
- Interaction dynamics
- Spatial and temporal analysis

**Usage:**
```python
from metainformant.simulation import agents

# Create simulation environment
world = agents.GridWorld(width=50, height=50, num_agents=100)

# Define agent behaviors
for agent in world.agents:
    agent.set_behavior(lambda: random.choice(["move", "interact"]))

# Run simulation
for _ in range(100):
    world.step()

# Analyze results
interactions = world.get_interaction_history()
```

## Advanced Simulation Features

### Evolutionary Dynamics
```python
from metainformant.simulation import sequences
from metainformant.math import price_equation

# Evolve sequences under selection
population = [sequences.generate_random_dna(100) for _ in range(100)]
for generation in range(50):
    # Calculate fitness
    fitness = [len(seq) for seq in population]  # Simple fitness function

    # Apply selection using Price equation
    selected = price_equation.apply_selection(population, fitness)

    # Generate next generation
    population = [sequences.mutate_sequence(seq) for seq in selected]
```

### Complex System Modeling
```python
from metainformant.simulation import agents, rna

# Multi-scale simulation
class GeneExpressionAgent(agents.Agent):
    def __init__(self, position):
        super().__init__(position)
        self.expression_level = 0

    def step(self):
        # Gene expression dynamics
        self.expression_level += rna.simulate_transcription_rate()
        self.move_towards_higher_expression()

# Run integrated simulation
world = agents.GridWorld(agents=[GeneExpressionAgent((i, j)) for i in range(10) for j in range(10)])
```

## Integration with Other Modules

### With Machine Learning Module
```python
from metainformant.simulation import sequences, rna
from metainformant.ml import classification

# Generate synthetic training data
sequences = [sequences.generate_random_dna(100) for _ in range(1000)]
labels = [0 if len(seq) < 500 else 1 for seq in sequences]

# Train classifier
model = classification.train_classifier(sequences, labels)
accuracy = classification.evaluate_model(model, sequences, labels)
```

### With Visualization Module
```python
from metainformant.simulation import agents
from metainformant.visualization import animations

# Visualize agent dynamics
world = agents.GridWorld(width=20, height=20, num_agents=50)
fig, anim = animations.animate_agent_movement(world, steps=100)
animations.save_animation(anim, "agent_simulation.gif")
```

## Performance Considerations

- **Memory Efficiency**: Streaming generation for large datasets
- **Parallel Generation**: Multi-threaded sequence and data generation
- **Configurable Complexity**: Scale from simple to complex simulations
- **Reproducibility**: Deterministic random seeds for reproducible results

## Parameter Control

All simulations support extensive parameter control:

```python
from metainformant.simulation import sequences

# Highly configurable sequence generation
params = {
    'length_distribution': 'normal',
    'length_mean': 1000,
    'length_std': 100,
    'gc_content': 0.5,
    'mutation_model': 'Jukes-Cantor'
}

sequences = sequences.generate_population(1000, **params)
```

## Testing and Validation

- **Distribution Validation**: Statistical tests for generated data distributions
- **Biological Realism**: Validation against real biological data
- **Performance Testing**: Benchmarking for large-scale simulations
- **Edge Case Handling**: Robust behavior with extreme parameters

## Dependencies

- **Core**: NumPy for numerical computations
- **Optional**: NetworkX for complex network simulations
- **Visualization**: Matplotlib for simulation visualization

## Usage Examples

### Synthetic Sequence Dataset
```python
from metainformant.simulation import sequences

# Generate comprehensive sequence dataset
dataset = sequences.generate_synthetic_dataset(
    n_sequences=10000,
    sequence_type='dna',
    length_range=(500, 2000),
    include_mutations=True,
    mutation_rate=0.001
)

# Save for downstream analysis
sequences.save_fasta(dataset['sequences'], "synthetic_dna.fasta")
```

### Agent-Based Ecosystem Simulation
```python
from metainformant.simulation import agents

# Ecosystem simulation
ecosystem = agents.Ecosystem(
    environment=agents.GridWorld(100, 100),
    species=[
        agents.PreySpecies(count=200, behavior="evasive"),
        agents.PredatorSpecies(count=20, behavior="pursuit")
    ]
)

# Run long-term simulation
for year in range(10):
    for day in range(365):
        ecosystem.step()
    ecosystem.analyze_annual_dynamics()
```

This module provides powerful tools for synthetic data generation and hypothesis testing, enabling researchers to explore biological systems in controlled computational environments.
