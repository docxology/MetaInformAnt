"""Agent-based ecosystem simulation utilities.

This module provides classes and functions for simulating ecosystems with
interacting agents, population dynamics, and evolutionary processes.
All simulations support reproducible results through random seed control.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import errors
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


@dataclass
class Agent:
    """Represents an individual agent in the ecosystem.

    Used by the Ecosystem-based simulation system for ecological modeling
    with typed agents (producer, consumer, decomposer).
    """

    id: int
    agent_type: str
    position: Tuple[int, int]
    properties: Dict[str, Any] = field(default_factory=dict)
    energy: float = 100.0
    age: int = 0
    alive: bool = True

    def move(self, new_position: Tuple[int, int]) -> None:
        """Move agent to new position."""
        self.position = new_position

    def update_energy(self, delta: float) -> None:
        """Update agent energy level."""
        self.energy += delta
        self.energy = max(0, self.energy)

    def reproduce(
        self, partner: Agent | None = None, mutation_rate: float = 0.01, rng: random.Random | None = None
    ) -> Agent:
        """Create offspring agent."""
        if rng is None:
            rng = random.Random()

        # Simple reproduction with inheritance and mutation
        offspring_properties = self.properties.copy()

        # Add mutations
        for key, value in offspring_properties.items():
            if isinstance(value, (int, float)) and rng.random() < mutation_rate:
                # Add small random variation
                variation = rng.normalvariate(0, abs(value) * 0.1) if value != 0 else rng.normalvariate(0, 1)
                offspring_properties[key] = value + variation

        return Agent(
            id=rng.randint(10000, 99999),  # Generate unique ID
            agent_type=self.agent_type,
            position=self.position,
            properties=offspring_properties,
            energy=self.energy * 0.5,  # Split energy with parent
        )

    def interact(self, other: Agent, interaction_type: str) -> Dict[str, Any]:
        """Interact with another agent."""
        # Base interaction - can be overridden by subclasses
        if interaction_type == "compete":
            # Competition reduces both agents' energy
            energy_loss = min(self.energy, other.energy) * 0.1
            self.update_energy(-energy_loss)
            other.update_energy(-energy_loss)
            return {"interaction": "competition", "energy_loss": energy_loss}

        elif interaction_type == "cooperate":
            # Cooperation benefits both
            energy_gain = (self.energy + other.energy) * 0.05
            self.update_energy(energy_gain)
            other.update_energy(energy_gain)
            return {"interaction": "cooperation", "energy_gain": energy_gain}

        return {"interaction": "neutral"}


@dataclass
class GridAgent:
    """Lightweight agent for grid-based simulations.

    Used by GridWorld for simple spatial agent-based models where agents
    have x/y coordinates and energy.
    """

    id: int
    x: int
    y: int
    energy: float = 100.0
    properties: Dict[str, Any] = field(default_factory=dict)

    def step(self, world: GridWorld, rng: random.Random) -> None:
        """Execute one movement step within a GridWorld.

        Args:
            world: The GridWorld environment
            rng: Random number generator
        """
        dx = rng.choice([-1, 0, 1])
        dy = rng.choice([-1, 0, 1])
        self.x = max(0, min(world.width - 1, self.x + dx))
        self.y = max(0, min(world.height - 1, self.y + dy))
        self.energy = max(0, self.energy - 0.1)


@dataclass
class Ecosystem:
    """Represents an ecosystem containing multiple agents."""

    size: Tuple[int, int]
    agents: Dict[int, Agent] = field(default_factory=dict)
    environment: Dict[str, Any] = field(default_factory=dict)
    time_step: int = 0

    def add_agent(self, agent: Agent) -> None:
        """Add an agent to the ecosystem."""
        self.agents[agent.id] = agent

    def remove_agent(self, agent_id: int) -> None:
        """Remove an agent from the ecosystem."""
        if agent_id in self.agents:
            del self.agents[agent_id]

    def get_agents_at_position(self, position: Tuple[int, int]) -> List[Agent]:
        """Get all agents at a specific position."""
        return [agent for agent in self.agents.values() if agent.position == position and agent.alive]

    def get_agents_by_type(self, agent_type: str) -> List[Agent]:
        """Get all agents of a specific type."""
        return [agent for agent in self.agents.values() if agent.agent_type == agent_type and agent.alive]

    def get_neighboring_positions(self, position: Tuple[int, int], radius: int = 1) -> List[Tuple[int, int]]:
        """Get positions within radius of given position."""
        x, y = position
        neighbors = []

        for dx in range(-radius, radius + 1):
            for dy in range(-radius, radius + 1):
                if dx == 0 and dy == 0:
                    continue
                nx, ny = x + dx, y + dy
                if 0 <= nx < self.size[0] and 0 <= ny < self.size[1]:
                    neighbors.append((nx, ny))

        return neighbors

    def get_random_position(self, rng: random.Random | None = None) -> Tuple[int, int]:
        """Get a random valid position in the ecosystem."""
        if rng is None:
            rng = random.Random()
        return (rng.randint(0, self.size[0] - 1), rng.randint(0, self.size[1] - 1))


def create_ecosystem(
    n_agents: int, agent_types: List[str], environment_size: Tuple[int, int], *, rng: random.Random | None = None
) -> Ecosystem:
    """Create an ecosystem with randomly placed agents.

    Args:
        n_agents: Total number of agents to create
        agent_types: List of possible agent types
        environment_size: Size of the environment (width, height)
        rng: Random number generator for reproducibility

    Returns:
        Initialized Ecosystem

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_agents, min_val=1, name="n_agents")
    validation.validate_type(agent_types, list, "agent_types")
    validation.validate_range(len(agent_types), min_val=1, name="agent_types length")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    ecosystem = Ecosystem(size=environment_size)

    # Create agents
    agents_per_type = n_agents // len(agent_types)
    remainder = n_agents % len(agent_types)

    agent_id = 1
    for i, agent_type in enumerate(agent_types):
        # Distribute agents evenly across types
        n_this_type = agents_per_type + (1 if i < remainder else 0)

        for _ in range(n_this_type):
            position = ecosystem.get_random_position(rng)

            # Create agent with type-specific properties
            properties = {}
            if agent_type == "producer":
                properties = {"reproduction_rate": 0.8, "energy_efficiency": 0.9}
            elif agent_type == "consumer":
                properties = {"hunting_skill": rng.uniform(0.5, 1.0), "energy_efficiency": 0.7}
            elif agent_type == "decomposer":
                properties = {"decomposition_rate": rng.uniform(0.6, 0.9), "energy_efficiency": 0.8}

            agent = Agent(
                id=agent_id,
                agent_type=agent_type,
                position=position,
                properties=properties,
                energy=rng.uniform(50, 150),
            )

            ecosystem.add_agent(agent)
            agent_id += 1

    return ecosystem


def run_simulation(
    ecosystem: Ecosystem, n_steps: int, *, record_interval: int = 1, rng: random.Random | None = None
) -> List[Dict[str, Any]]:
    """Run ecosystem simulation for specified number of steps.

    Args:
        ecosystem: The ecosystem to simulate
        n_steps: Number of simulation steps
        record_interval: Interval at which to record state
        rng: Random number generator for reproducibility

    Returns:
        List of simulation snapshots

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_steps, min_val=1, name="n_steps")
    validation.validate_range(record_interval, min_val=1, name="record_interval")

    if rng is None:
        rng = random.Random()

    snapshots = []
    logger.info(f"Starting ecosystem simulation for {n_steps} steps")

    for step in range(n_steps):
        # Record state at intervals
        if step % record_interval == 0:
            snapshot = {
                "step": step,
                "n_agents": len([a for a in ecosystem.agents.values() if a.alive]),
                "agent_counts": count_agents_by_type(ecosystem),
                "total_energy": sum(a.energy for a in ecosystem.agents.values() if a.alive),
                "agent_positions": [(a.id, a.position) for a in ecosystem.agents.values() if a.alive],
            }
            snapshots.append(snapshot)

        # Run one simulation step
        simulation_step(ecosystem, rng)

        ecosystem.time_step += 1

    logger.info(f"Simulation completed. Recorded {len(snapshots)} snapshots")
    return snapshots


def simulation_step(ecosystem: Ecosystem, rng: random.Random) -> None:
    """Execute one step of the ecosystem simulation."""
    # Process each agent
    agents_to_process = list(ecosystem.agents.values())
    rng.shuffle(agents_to_process)  # Randomize order

    new_agents = []

    for agent in agents_to_process:
        if not agent.alive:
            continue

        # Age the agent
        agent.age += 1

        # Energy consumption (basal metabolism)
        agent.update_energy(-1.0)

        # Agent behavior based on type
        if agent.agent_type == "producer":
            # Producers gain energy from environment
            energy_gain = 2.0 * agent.properties.get("energy_efficiency", 1.0)
            agent.update_energy(energy_gain)

        elif agent.agent_type == "consumer":
            # Consumers hunt for food
            neighbors = ecosystem.get_neighboring_positions(agent.position, radius=2)
            prey_found = False

            for pos in neighbors:
                nearby_agents = ecosystem.get_agents_at_position(pos)
                for nearby in nearby_agents:
                    if nearby.agent_type == "producer" and nearby.alive:
                        # Successful hunt
                        hunting_skill = agent.properties.get("hunting_skill", 0.5)
                        if rng.random() < hunting_skill:
                            energy_gain = nearby.energy * 0.5
                            agent.update_energy(energy_gain)
                            nearby.alive = False  # Prey dies
                            prey_found = True
                            break
                if prey_found:
                    break

        elif agent.agent_type == "decomposer":
            # Decomposers consume dead matter
            neighbors = ecosystem.get_neighboring_positions(agent.position, radius=1)
            dead_found = False

            for pos in neighbors:
                nearby_agents = ecosystem.get_agents_at_position(pos)
                for nearby in nearby_agents:
                    if not nearby.alive:
                        # Consume dead agent
                        decomposition_rate = agent.properties.get("decomposition_rate", 0.7)
                        energy_gain = nearby.energy * decomposition_rate
                        agent.update_energy(energy_gain)
                        ecosystem.remove_agent(nearby.id)  # Remove corpse
                        dead_found = True
                        break
                if dead_found:
                    break

        # Movement
        if rng.random() < 0.3:  # 30% chance to move each step
            new_position = ecosystem.get_random_position(rng)
            agent.move(new_position)

        # Reproduction
        reproduction_chance = 0.05  # Base reproduction chance
        if agent.agent_type == "producer":
            reproduction_chance *= agent.properties.get("reproduction_rate", 1.0)

        if agent.energy > 150 and rng.random() < reproduction_chance:
            offspring = agent.reproduce(mutation_rate=0.01, rng=rng)
            offspring.move(ecosystem.get_random_position(rng))
            new_agents.append(offspring)
            agent.update_energy(-50)  # Energy cost of reproduction

        # Death from starvation
        if agent.energy <= 0:
            agent.alive = False
            continue

        # Age-related death
        if agent.age > 100 and rng.random() < 0.01:
            agent.alive = False
            continue

    # Add new agents to ecosystem
    for new_agent in new_agents:
        ecosystem.add_agent(new_agent)


def add_agent(
    ecosystem: Ecosystem, agent_type: str, position: Tuple[int, int], *, properties: Dict[str, Any] | None = None
) -> None:
    """Add a new agent to the ecosystem.

    Args:
        ecosystem: The ecosystem to modify
        agent_type: Type of agent to add
        position: Position to place the agent
        properties: Optional agent properties

    Raises:
        ValueError: If position is invalid
    """
    x, y = position
    if not (0 <= x < ecosystem.size[0] and 0 <= y < ecosystem.size[1]):
        raise errors.ValidationError(f"Position {position} is outside ecosystem bounds {ecosystem.size}")

    if properties is None:
        properties = {}

    # Generate unique ID
    agent_id = max(ecosystem.agents.keys()) + 1 if ecosystem.agents else 1

    agent = Agent(id=agent_id, agent_type=agent_type, position=position, properties=properties, energy=100.0)

    ecosystem.add_agent(agent)


def remove_agent(ecosystem: Ecosystem, agent_id: int) -> None:
    """Remove an agent from the ecosystem.

    Args:
        ecosystem: The ecosystem to modify
        agent_id: ID of agent to remove
    """
    ecosystem.remove_agent(agent_id)


def get_population_dynamics(simulation_data: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Analyze population dynamics from simulation snapshots.

    Args:
        simulation_data: List of simulation snapshots

    Returns:
        Dictionary with population dynamics analysis
    """
    if not simulation_data:
        return {}

    steps = [d["step"] for d in simulation_data]
    total_agents = [d["n_agents"] for d in simulation_data]

    # Analyze agent type dynamics
    agent_types = set()
    for data in simulation_data:
        agent_types.update(data["agent_counts"].keys())

    type_dynamics = {}
    for agent_type in agent_types:
        counts = [d["agent_counts"].get(agent_type, 0) for d in simulation_data]
        type_dynamics[agent_type] = {
            "counts": counts,
            "max_count": max(counts),
            "min_count": min(counts),
            "final_count": counts[-1],
        }

    # Calculate diversity metrics
    diversity_indices = []
    for data in simulation_data:
        counts = list(data["agent_counts"].values())
        if counts:
            # Shannon diversity index
            total = sum(counts)
            proportions = [c / total for c in counts]
            shannon = -sum(p * np.log(p) for p in proportions if p > 0)
            diversity_indices.append(shannon)
        else:
            diversity_indices.append(0)

    return {
        "total_agents_trajectory": total_agents,
        "agent_type_dynamics": type_dynamics,
        "diversity_trajectory": diversity_indices,
        "simulation_steps": steps,
        "final_diversity": diversity_indices[-1] if diversity_indices else 0,
        "max_diversity": max(diversity_indices) if diversity_indices else 0,
    }


def calculate_biodiversity_metrics(simulation_data: List[Dict[str, Any]]) -> Dict[str, float]:
    """Calculate biodiversity metrics from simulation data.

    Args:
        simulation_data: List of simulation snapshots

    Returns:
        Dictionary with biodiversity metrics
    """
    if not simulation_data:
        return {}

    final_snapshot = simulation_data[-1]
    agent_counts = final_snapshot["agent_counts"]

    total_individuals = sum(agent_counts.values())
    n_types = len(agent_counts)

    # Species richness
    richness = n_types

    # Shannon diversity
    if total_individuals > 0:
        proportions = [count / total_individuals for count in agent_counts.values()]
        shannon = -sum(p * np.log(p) for p in proportions if p > 0)
    else:
        shannon = 0

    # Simpson diversity
    if total_individuals > 0:
        proportions = [count / total_individuals for count in agent_counts.values()]
        simpson = 1 - sum(p**2 for p in proportions)
    else:
        simpson = 0

    # Evenness (Pielou's evenness)
    if n_types > 1 and shannon > 0:
        evenness = shannon / np.log(n_types)
    else:
        evenness = 0

    return {
        "species_richness": richness,
        "shannon_diversity": shannon,
        "simpson_diversity": simpson,
        "pielou_evenness": evenness,
        "total_individuals": total_individuals,
    }


def count_agents_by_type(ecosystem: Ecosystem) -> Dict[str, int]:
    """Count living agents by type.

    Args:
        ecosystem: The ecosystem to analyze

    Returns:
        Dictionary mapping agent types to counts
    """
    counts = {}
    for agent in ecosystem.agents.values():
        if agent.alive:
            counts[agent.agent_type] = counts.get(agent.agent_type, 0) + 1
    return counts


def simulate_predator_prey(
    ecosystem: Ecosystem,
    n_steps: int,
    *,
    predator_efficiency: float = 0.8,
    prey_reproduction: float = 0.1,
    rng: random.Random | None = None,
) -> List[Dict[str, Any]]:
    """Run a predator-prey simulation.

    Args:
        ecosystem: Initialized ecosystem
        n_steps: Number of simulation steps
        predator_efficiency: Efficiency of predation (0-1)
        prey_reproduction: Prey reproduction rate
        rng: Random number generator for reproducibility

    Returns:
        List of simulation snapshots with predator-prey dynamics
    """
    validation.validate_range(predator_efficiency, min_val=0.0, max_val=1.0, name="predator_efficiency")
    validation.validate_range(prey_reproduction, min_val=0.0, name="prey_reproduction")

    if rng is None:
        rng = random.Random()

    snapshots = []
    logger.info(f"Starting predator-prey simulation for {n_steps} steps")

    for step in range(n_steps):
        # Custom predator-prey dynamics
        predators = ecosystem.get_agents_by_type("consumer")
        prey = ecosystem.get_agents_by_type("producer")

        # Predators hunt prey
        for predator in predators:
            if not predator.alive:
                continue

            # Find nearby prey
            neighbors = ecosystem.get_neighboring_positions(predator.position, radius=3)
            nearby_prey = []

            for pos in neighbors:
                nearby_prey.extend(ecosystem.get_agents_at_position(pos))

            nearby_prey = [a for a in nearby_prey if a.agent_type == "producer" and a.alive]

            if nearby_prey and rng.random() < predator_efficiency:
                # Successful hunt
                target_prey = rng.choice(nearby_prey)
                energy_gain = target_prey.energy * 0.7
                predator.update_energy(energy_gain)
                target_prey.alive = False

        # Prey reproduction
        for prey_agent in prey:
            if prey_agent.alive and prey_agent.energy > 80 and rng.random() < prey_reproduction:
                # Create offspring
                offspring = prey_agent.reproduce(rng=rng)
                offspring.move(ecosystem.get_random_position(rng))
                ecosystem.add_agent(offspring)
                prey_agent.update_energy(-30)

        # Standard simulation step
        simulation_step(ecosystem, rng)

        # Record state
        snapshot = {
            "step": step,
            "n_predators": len(predators),
            "n_prey": len(prey),
            "total_energy": sum(a.energy for a in ecosystem.agents.values() if a.alive),
        }
        snapshots.append(snapshot)

        ecosystem.time_step += 1

    logger.info("Predator-prey simulation completed")
    return snapshots


def simulate_competition(
    ecosystem: Ecosystem, n_steps: int, *, competition_strength: float = 0.5, rng: random.Random | None = None
) -> List[Dict[str, Any]]:
    """Run a competition simulation between agent types.

    Args:
        ecosystem: Initialized ecosystem
        n_steps: Number of simulation steps
        competition_strength: Strength of competitive interactions
        rng: Random number generator for reproducibility

    Returns:
        List of simulation snapshots
    """
    validation.validate_range(competition_strength, min_val=0.0, max_val=1.0, name="competition_strength")

    if rng is None:
        rng = random.Random()

    snapshots = []
    logger.info(f"Starting competition simulation for {n_steps} steps")

    for step in range(n_steps):
        # Competition between agents at same position
        positions_occupied = set()

        for agent in list(ecosystem.agents.values()):
            if not agent.alive:
                continue

            positions_occupied.add(agent.position)

            # Check for competitors at same position
            competitors = ecosystem.get_agents_at_position(agent.position)
            competitors = [c for c in competitors if c.id != agent.id and c.alive]

            for competitor in competitors:
                if rng.random() < competition_strength:
                    # Competition reduces energy of both
                    energy_loss = min(agent.energy, competitor.energy) * 0.2
                    agent.update_energy(-energy_loss)
                    competitor.update_energy(-energy_loss)

        # Standard simulation step
        simulation_step(ecosystem, rng)

        # Record state
        snapshot = {
            "step": step,
            "occupied_positions": len(positions_occupied),
            "agent_counts": count_agents_by_type(ecosystem),
            "competition_events": len(
                [a for a in ecosystem.agents.values() if not a.alive and a.energy <= 0]
            ),  # Deaths from competition
        }
        snapshots.append(snapshot)

        ecosystem.time_step += 1

    logger.info("Competition simulation completed")
    return snapshots


class GridWorld:
    """A grid-based world for agent-based simulations.

    This class provides a 2D grid environment where GridAgents can move,
    interact, and evolve according to specified rules.

    Args:
        width: Width of the grid
        height: Height of the grid
        num_agents: Number of agents to create
        rng: Random number generator (optional)
        **kwargs: Additional parameters
    """

    def __init__(self, width: int, height: int, num_agents: int, rng: Optional[Any] = None, **kwargs: Any):
        self.width = width
        self.height = height
        self.agents: List[GridAgent] = []

        if rng is None:
            self.rng = random.Random()
        else:
            self.rng = rng

        for i in range(num_agents):
            agent = GridAgent(
                id=i,
                x=self.rng.randint(0, width - 1),
                y=self.rng.randint(0, height - 1),
                energy=100.0,
            )
            self.agents.append(agent)

        self.time_step = 0

    def step(self) -> Dict[str, Any]:
        """Execute one simulation step.

        Returns:
            Dictionary with step results
        """
        for agent in self.agents:
            agent.step(self, self.rng)

        self.time_step += 1

        return {
            "time_step": self.time_step,
            "num_agents": len(self.agents),
            "agent_positions": [(a.x, a.y) for a in self.agents],
            "total_energy": sum(a.energy for a in self.agents),
        }

    def positions(self) -> List[Tuple[int, int]]:
        """Get current positions of all agents.

        Returns:
            List of (x, y) coordinate tuples
        """
        return [(agent.x, agent.y) for agent in self.agents]

    def get_agent_positions(self) -> List[Tuple[int, int]]:
        """Get current positions of all agents (alias for positions()).

        Returns:
            List of (x, y) coordinate tuples
        """
        return self.positions()

    def get_agent_energies(self) -> List[float]:
        """Get current energy levels of all agents.

        Returns:
            List of energy values
        """
        return [agent.energy for agent in self.agents]

    def add_agent(self, x: Optional[int] = None, y: Optional[int] = None, energy: float = 100.0) -> GridAgent:
        """Add a new agent to the world.

        Args:
            x: X coordinate (random if None)
            y: Y coordinate (random if None)
            energy: Initial energy level

        Returns:
            The newly created GridAgent
        """
        if x is None:
            x = self.rng.randint(0, self.width - 1)
        if y is None:
            y = self.rng.randint(0, self.height - 1)

        agent = GridAgent(id=len(self.agents), x=x, y=y, energy=energy)
        self.agents.append(agent)
        return agent

    def remove_agent(self, agent: GridAgent) -> None:
        """Remove an agent from the world.

        Args:
            agent: GridAgent to remove
        """
        if agent in self.agents:
            self.agents.remove(agent)

    def get_neighbors(self, agent: GridAgent, radius: int = 1) -> List[GridAgent]:
        """Get neighboring agents within a given radius (Manhattan distance).

        Args:
            agent: Reference agent
            radius: Search radius

        Returns:
            List of neighboring GridAgents
        """
        neighbors = []
        for other in self.agents:
            if other is agent:
                continue
            distance = abs(other.x - agent.x) + abs(other.y - agent.y)
            if distance <= radius:
                neighbors.append(other)
        return neighbors

    def run_simulation(self, num_steps: int) -> List[Dict[str, Any]]:
        """Run simulation for multiple steps.

        Args:
            num_steps: Number of steps to run

        Returns:
            List of step results
        """
        results = []
        for _ in range(num_steps):
            step_result = self.step()
            results.append(step_result)
        return results
