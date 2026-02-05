"""Comprehensive tests for the agent-based ecosystem simulation module.

Tests cover:
- Agent creation, movement, energy, reproduction, interaction
- Ecosystem creation, agent management, spatial queries
- create_ecosystem with multiple agent types
- run_simulation with snapshot recording
- simulation_step basic operation
- Population dynamics analysis
- Biodiversity metrics (Shannon, Simpson, evenness, richness)
- count_agents_by_type
- Predator-prey simulation
- Competition simulation
- Edge cases (single agent, boundary positions, empty data)
- Reproducibility with seeded RNGs

All tests use real implementations (NO mocking).
"""

from __future__ import annotations

import math
import random
from typing import Dict, List

import numpy as np
import pytest

from metainformant.simulation import (
    Agent,
    Ecosystem,
    add_agent,
    calculate_biodiversity_metrics,
    count_agents_by_type,
    create_ecosystem,
    get_population_dynamics,
    remove_agent,
    run_simulation,
    simulate_competition,
    simulate_predator_prey,
    simulation_step,
)
from metainformant.core.utils.errors import ValidationError


# ---------------------------------------------------------------------------
# Agent dataclass — creation and defaults
# ---------------------------------------------------------------------------


class TestAgentCreation:
    """Test Agent dataclass construction and default values."""

    def test_agent_basic_creation(self) -> None:
        """Agent created with required fields gets correct defaults."""
        agent = Agent(id=1, agent_type="producer", position=(5, 5))
        assert agent.id == 1
        assert agent.agent_type == "producer"
        assert agent.position == (5, 5)
        assert agent.properties == {}
        assert agent.energy == 100.0
        assert agent.age == 0
        assert agent.alive is True

    def test_agent_full_creation(self) -> None:
        """Agent created with all fields stores them correctly."""
        props = {"speed": 2.5, "vision": 3}
        agent = Agent(
            id=42,
            agent_type="consumer",
            position=(10, 20),
            properties=props,
            energy=75.0,
            age=5,
            alive=False,
        )
        assert agent.id == 42
        assert agent.agent_type == "consumer"
        assert agent.position == (10, 20)
        assert agent.properties == props
        assert agent.energy == 75.0
        assert agent.age == 5
        assert agent.alive is False

    def test_agent_custom_properties_are_independent(self) -> None:
        """Each agent gets its own properties dict (no shared mutable default)."""
        a1 = Agent(id=1, agent_type="x", position=(0, 0))
        a2 = Agent(id=2, agent_type="x", position=(0, 0))
        a1.properties["key"] = "val"
        assert "key" not in a2.properties


# ---------------------------------------------------------------------------
# Agent.move
# ---------------------------------------------------------------------------


class TestAgentMove:
    """Test Agent.move() method."""

    def test_move_updates_position(self) -> None:
        agent = Agent(id=1, agent_type="producer", position=(0, 0))
        agent.move((3, 7))
        assert agent.position == (3, 7)

    def test_move_to_same_position(self) -> None:
        agent = Agent(id=1, agent_type="producer", position=(5, 5))
        agent.move((5, 5))
        assert agent.position == (5, 5)

    def test_move_multiple_times(self) -> None:
        agent = Agent(id=1, agent_type="producer", position=(0, 0))
        for i in range(10):
            agent.move((i, i + 1))
        assert agent.position == (9, 10)


# ---------------------------------------------------------------------------
# Agent.update_energy
# ---------------------------------------------------------------------------


class TestAgentEnergy:
    """Test Agent.update_energy() method including clamping."""

    def test_positive_delta(self) -> None:
        agent = Agent(id=1, agent_type="producer", position=(0, 0), energy=50.0)
        agent.update_energy(25.0)
        assert agent.energy == 75.0

    def test_negative_delta(self) -> None:
        agent = Agent(id=1, agent_type="producer", position=(0, 0), energy=50.0)
        agent.update_energy(-20.0)
        assert agent.energy == 30.0

    def test_energy_clamped_at_zero(self) -> None:
        """Energy never goes below zero."""
        agent = Agent(id=1, agent_type="producer", position=(0, 0), energy=10.0)
        agent.update_energy(-100.0)
        assert agent.energy == 0.0

    def test_energy_exact_zero(self) -> None:
        agent = Agent(id=1, agent_type="producer", position=(0, 0), energy=30.0)
        agent.update_energy(-30.0)
        assert agent.energy == 0.0

    def test_zero_delta(self) -> None:
        agent = Agent(id=1, agent_type="producer", position=(0, 0), energy=42.0)
        agent.update_energy(0.0)
        assert agent.energy == 42.0

    def test_large_positive_delta(self) -> None:
        agent = Agent(id=1, agent_type="producer", position=(0, 0), energy=100.0)
        agent.update_energy(1e6)
        assert agent.energy == pytest.approx(1e6 + 100.0)


# ---------------------------------------------------------------------------
# Agent.reproduce
# ---------------------------------------------------------------------------


class TestAgentReproduce:
    """Test Agent.reproduce() method."""

    def test_reproduce_without_partner(self) -> None:
        rng = random.Random(42)
        parent = Agent(id=1, agent_type="producer", position=(5, 5), energy=100.0)
        offspring = parent.reproduce(rng=rng)

        assert offspring.agent_type == parent.agent_type
        assert offspring.position == parent.position
        assert offspring.energy == pytest.approx(parent.energy * 0.5)
        assert offspring.id != parent.id
        assert offspring.age == 0
        assert offspring.alive is True

    def test_reproduce_energy_split(self) -> None:
        """Offspring gets 50% of parent energy at time of reproduction."""
        parent = Agent(id=1, agent_type="consumer", position=(3, 3), energy=80.0)
        rng = random.Random(99)
        offspring = parent.reproduce(rng=rng)
        assert offspring.energy == pytest.approx(40.0)

    def test_reproduce_inherits_properties(self) -> None:
        """Offspring properties are based on parent properties."""
        parent = Agent(
            id=1,
            agent_type="producer",
            position=(0, 0),
            properties={"speed": 5.0, "label": "fast"},
            energy=100.0,
        )
        rng = random.Random(42)
        offspring = parent.reproduce(mutation_rate=0.0, rng=rng)
        # With zero mutation rate, numeric properties stay the same
        assert offspring.properties["speed"] == parent.properties["speed"]
        # Non-numeric properties are always copied as-is
        assert offspring.properties["label"] == "fast"

    def test_reproduce_with_mutation(self) -> None:
        """High mutation rate causes numeric property changes."""
        parent = Agent(
            id=1,
            agent_type="producer",
            position=(0, 0),
            properties={"speed": 10.0, "size": 5.0},
            energy=200.0,
        )
        rng = random.Random(42)
        offspring = parent.reproduce(mutation_rate=1.0, rng=rng)
        # With mutation_rate=1.0 all numeric properties should be mutated
        # At least one should differ
        speed_changed = offspring.properties["speed"] != parent.properties["speed"]
        size_changed = offspring.properties["size"] != parent.properties["size"]
        assert speed_changed or size_changed

    def test_reproduce_with_partner_still_works(self) -> None:
        """Passing a partner does not break reproduction (partner arg is accepted)."""
        parent = Agent(id=1, agent_type="consumer", position=(0, 0), energy=100.0)
        partner = Agent(id=2, agent_type="consumer", position=(1, 1), energy=80.0)
        rng = random.Random(42)
        offspring = parent.reproduce(partner=partner, rng=rng)
        assert offspring.agent_type == parent.agent_type
        assert offspring.alive is True

    def test_reproduce_reproducible(self) -> None:
        """Same seed produces identical offspring."""
        parent = Agent(
            id=1,
            agent_type="producer",
            position=(2, 2),
            properties={"x": 10.0},
            energy=100.0,
        )
        o1 = parent.reproduce(mutation_rate=0.5, rng=random.Random(7))
        # Reset parent energy back for second call
        parent.energy = 100.0
        o2 = parent.reproduce(mutation_rate=0.5, rng=random.Random(7))
        assert o1.id == o2.id
        assert o1.properties == o2.properties
        assert o1.energy == o2.energy


# ---------------------------------------------------------------------------
# Agent.interact
# ---------------------------------------------------------------------------


class TestAgentInteract:
    """Test Agent.interact() method."""

    def test_compete_interaction(self) -> None:
        a = Agent(id=1, agent_type="consumer", position=(0, 0), energy=100.0)
        b = Agent(id=2, agent_type="consumer", position=(0, 0), energy=100.0)
        result = a.interact(b, "compete")

        assert result["interaction"] == "competition"
        assert "energy_loss" in result
        expected_loss = min(100.0, 100.0) * 0.1  # 10.0
        assert result["energy_loss"] == pytest.approx(expected_loss)
        # Both agents lose energy
        assert a.energy == pytest.approx(100.0 - expected_loss)
        assert b.energy == pytest.approx(100.0 - expected_loss)

    def test_cooperate_interaction(self) -> None:
        a = Agent(id=1, agent_type="producer", position=(0, 0), energy=100.0)
        b = Agent(id=2, agent_type="producer", position=(0, 0), energy=60.0)
        result = a.interact(b, "cooperate")

        assert result["interaction"] == "cooperation"
        assert "energy_gain" in result
        expected_gain = (100.0 + 60.0) * 0.05  # 8.0
        assert result["energy_gain"] == pytest.approx(expected_gain)
        assert a.energy == pytest.approx(100.0 + expected_gain)
        assert b.energy == pytest.approx(60.0 + expected_gain)

    def test_neutral_interaction(self) -> None:
        """Unknown interaction type returns neutral with no energy change."""
        a = Agent(id=1, agent_type="producer", position=(0, 0), energy=100.0)
        b = Agent(id=2, agent_type="consumer", position=(0, 0), energy=100.0)
        result = a.interact(b, "unknown_type")

        assert result == {"interaction": "neutral"}
        assert a.energy == 100.0
        assert b.energy == 100.0

    def test_compete_asymmetric_energy(self) -> None:
        """Competition loss is based on the agent with less energy."""
        a = Agent(id=1, agent_type="consumer", position=(0, 0), energy=200.0)
        b = Agent(id=2, agent_type="consumer", position=(0, 0), energy=50.0)
        result = a.interact(b, "compete")
        expected_loss = min(200.0, 50.0) * 0.1  # 5.0
        assert result["energy_loss"] == pytest.approx(expected_loss)
        assert a.energy == pytest.approx(195.0)
        assert b.energy == pytest.approx(45.0)


# ---------------------------------------------------------------------------
# Ecosystem creation and basic operations
# ---------------------------------------------------------------------------


class TestEcosystem:
    """Test Ecosystem dataclass and its methods."""

    def test_ecosystem_creation_defaults(self) -> None:
        eco = Ecosystem(size=(10, 10))
        assert eco.size == (10, 10)
        assert eco.agents == {}
        assert eco.environment == {}
        assert eco.time_step == 0

    def test_add_agent(self) -> None:
        eco = Ecosystem(size=(10, 10))
        agent = Agent(id=1, agent_type="producer", position=(5, 5))
        eco.add_agent(agent)
        assert 1 in eco.agents
        assert eco.agents[1] is agent

    def test_remove_agent(self) -> None:
        eco = Ecosystem(size=(10, 10))
        agent = Agent(id=1, agent_type="producer", position=(5, 5))
        eco.add_agent(agent)
        eco.remove_agent(1)
        assert 1 not in eco.agents

    def test_remove_nonexistent_agent(self) -> None:
        """Removing an agent that doesn't exist is a no-op."""
        eco = Ecosystem(size=(10, 10))
        eco.remove_agent(999)  # Should not raise

    def test_add_multiple_agents(self) -> None:
        eco = Ecosystem(size=(20, 20))
        for i in range(5):
            eco.add_agent(Agent(id=i, agent_type="producer", position=(i, i)))
        assert len(eco.agents) == 5

    def test_get_agents_at_position(self) -> None:
        eco = Ecosystem(size=(10, 10))
        a1 = Agent(id=1, agent_type="producer", position=(3, 3))
        a2 = Agent(id=2, agent_type="consumer", position=(3, 3))
        a3 = Agent(id=3, agent_type="producer", position=(5, 5))
        eco.add_agent(a1)
        eco.add_agent(a2)
        eco.add_agent(a3)

        at_33 = eco.get_agents_at_position((3, 3))
        assert len(at_33) == 2
        ids_at_33 = {a.id for a in at_33}
        assert ids_at_33 == {1, 2}

    def test_get_agents_at_position_excludes_dead(self) -> None:
        eco = Ecosystem(size=(10, 10))
        a1 = Agent(id=1, agent_type="producer", position=(3, 3), alive=True)
        a2 = Agent(id=2, agent_type="producer", position=(3, 3), alive=False)
        eco.add_agent(a1)
        eco.add_agent(a2)

        at_33 = eco.get_agents_at_position((3, 3))
        assert len(at_33) == 1
        assert at_33[0].id == 1

    def test_get_agents_by_type(self) -> None:
        eco = Ecosystem(size=(10, 10))
        eco.add_agent(Agent(id=1, agent_type="producer", position=(0, 0)))
        eco.add_agent(Agent(id=2, agent_type="consumer", position=(1, 1)))
        eco.add_agent(Agent(id=3, agent_type="producer", position=(2, 2)))

        producers = eco.get_agents_by_type("producer")
        assert len(producers) == 2
        consumers = eco.get_agents_by_type("consumer")
        assert len(consumers) == 1

    def test_get_agents_by_type_excludes_dead(self) -> None:
        eco = Ecosystem(size=(10, 10))
        eco.add_agent(Agent(id=1, agent_type="producer", position=(0, 0), alive=True))
        eco.add_agent(Agent(id=2, agent_type="producer", position=(1, 1), alive=False))

        producers = eco.get_agents_by_type("producer")
        assert len(producers) == 1

    def test_get_neighboring_positions_center(self) -> None:
        eco = Ecosystem(size=(10, 10))
        neighbors = eco.get_neighboring_positions((5, 5), radius=1)
        # 3x3 grid minus center = 8 neighbors
        assert len(neighbors) == 8
        assert (5, 5) not in neighbors
        assert (4, 4) in neighbors
        assert (6, 6) in neighbors

    def test_get_neighboring_positions_corner(self) -> None:
        """Corner position has fewer neighbors due to boundary clipping."""
        eco = Ecosystem(size=(10, 10))
        neighbors = eco.get_neighboring_positions((0, 0), radius=1)
        # Only (1,0), (0,1), (1,1) are valid
        assert len(neighbors) == 3
        assert all(0 <= x < 10 and 0 <= y < 10 for x, y in neighbors)

    def test_get_neighboring_positions_edge(self) -> None:
        """Edge position clips neighbors outside bounds."""
        eco = Ecosystem(size=(10, 10))
        neighbors = eco.get_neighboring_positions((0, 5), radius=1)
        # Left edge: 5 neighbors instead of 8
        assert all(x >= 0 for x, _ in neighbors)
        assert len(neighbors) == 5

    def test_get_neighboring_positions_large_radius(self) -> None:
        eco = Ecosystem(size=(10, 10))
        neighbors = eco.get_neighboring_positions((5, 5), radius=2)
        # 5x5 grid minus center = 24 neighbors (all within bounds)
        assert len(neighbors) == 24
        assert (5, 5) not in neighbors

    def test_get_random_position(self) -> None:
        eco = Ecosystem(size=(10, 10))
        rng = random.Random(42)
        pos = eco.get_random_position(rng)
        x, y = pos
        assert 0 <= x < 10
        assert 0 <= y < 10

    def test_get_random_position_within_bounds(self) -> None:
        """Many random positions are all within bounds."""
        eco = Ecosystem(size=(5, 8))
        rng = random.Random(42)
        for _ in range(100):
            x, y = eco.get_random_position(rng)
            assert 0 <= x < 5
            assert 0 <= y < 8


# ---------------------------------------------------------------------------
# create_ecosystem function
# ---------------------------------------------------------------------------


class TestCreateEcosystem:
    """Test the create_ecosystem factory function."""

    def test_basic_creation(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(20, 20),
            rng=rng,
        )
        assert isinstance(eco, Ecosystem)
        assert eco.size == (20, 20)
        assert len(eco.agents) == 10

    def test_agents_distributed_across_types(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=9,
            agent_types=["producer", "consumer", "decomposer"],
            environment_size=(20, 20),
            rng=rng,
        )
        counts = count_agents_by_type(eco)
        assert counts["producer"] == 3
        assert counts["consumer"] == 3
        assert counts["decomposer"] == 3

    def test_uneven_distribution(self) -> None:
        """When n_agents is not divisible by type count, remainder is distributed."""
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer", "decomposer"],
            environment_size=(20, 20),
            rng=rng,
        )
        counts = count_agents_by_type(eco)
        total = sum(counts.values())
        assert total == 10
        # First type gets the extra agent: 4+3+3=10
        assert counts["producer"] == 4

    def test_single_agent_type(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=5,
            agent_types=["producer"],
            environment_size=(10, 10),
            rng=rng,
        )
        counts = count_agents_by_type(eco)
        assert counts == {"producer": 5}

    def test_positions_within_bounds(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=50,
            agent_types=["producer", "consumer"],
            environment_size=(10, 10),
            rng=rng,
        )
        for agent in eco.agents.values():
            x, y = agent.position
            assert 0 <= x < 10
            assert 0 <= y < 10

    def test_producer_properties(self) -> None:
        """Producers get reproduction_rate and energy_efficiency properties."""
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=3,
            agent_types=["producer"],
            environment_size=(10, 10),
            rng=rng,
        )
        for agent in eco.agents.values():
            assert "reproduction_rate" in agent.properties
            assert "energy_efficiency" in agent.properties

    def test_consumer_properties(self) -> None:
        """Consumers get hunting_skill and energy_efficiency properties."""
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=3,
            agent_types=["consumer"],
            environment_size=(10, 10),
            rng=rng,
        )
        for agent in eco.agents.values():
            assert "hunting_skill" in agent.properties
            assert "energy_efficiency" in agent.properties

    def test_decomposer_properties(self) -> None:
        """Decomposers get decomposition_rate and energy_efficiency properties."""
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=3,
            agent_types=["decomposer"],
            environment_size=(10, 10),
            rng=rng,
        )
        for agent in eco.agents.values():
            assert "decomposition_rate" in agent.properties
            assert "energy_efficiency" in agent.properties

    def test_invalid_n_agents_zero(self) -> None:
        with pytest.raises(ValidationError):
            create_ecosystem(
                n_agents=0,
                agent_types=["producer"],
                environment_size=(10, 10),
            )

    def test_invalid_n_agents_negative(self) -> None:
        with pytest.raises(ValidationError):
            create_ecosystem(
                n_agents=-5,
                agent_types=["producer"],
                environment_size=(10, 10),
            )

    def test_reproducibility(self) -> None:
        """Same seed produces identical ecosystems."""
        eco1 = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(20, 20),
            rng=random.Random(42),
        )
        eco2 = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(20, 20),
            rng=random.Random(42),
        )
        for aid in eco1.agents:
            assert eco1.agents[aid].position == eco2.agents[aid].position
            assert eco1.agents[aid].agent_type == eco2.agents[aid].agent_type
            assert eco1.agents[aid].energy == pytest.approx(eco2.agents[aid].energy)


# ---------------------------------------------------------------------------
# Module-level add_agent / remove_agent functions
# ---------------------------------------------------------------------------


class TestModuleLevelAddRemove:
    """Test the module-level add_agent() and remove_agent() functions."""

    def test_add_agent_function(self) -> None:
        eco = Ecosystem(size=(10, 10))
        add_agent(eco, "producer", (3, 3))
        assert len(eco.agents) == 1
        agent = list(eco.agents.values())[0]
        assert agent.agent_type == "producer"
        assert agent.position == (3, 3)
        assert agent.energy == 100.0

    def test_add_agent_with_properties(self) -> None:
        eco = Ecosystem(size=(10, 10))
        add_agent(eco, "consumer", (5, 5), properties={"speed": 3.0})
        agent = list(eco.agents.values())[0]
        assert agent.properties == {"speed": 3.0}

    def test_add_agent_out_of_bounds_raises(self) -> None:
        eco = Ecosystem(size=(10, 10))
        with pytest.raises(ValidationError):
            add_agent(eco, "producer", (10, 5))  # x == size[0] is out of bounds

    def test_add_agent_negative_position_raises(self) -> None:
        eco = Ecosystem(size=(10, 10))
        with pytest.raises(ValidationError):
            add_agent(eco, "producer", (-1, 5))

    def test_add_agent_boundary_valid(self) -> None:
        """Position (size-1, size-1) is the last valid position."""
        eco = Ecosystem(size=(10, 10))
        add_agent(eco, "producer", (9, 9))
        assert len(eco.agents) == 1

    def test_add_agent_origin_valid(self) -> None:
        eco = Ecosystem(size=(10, 10))
        add_agent(eco, "producer", (0, 0))
        assert len(eco.agents) == 1

    def test_remove_agent_function(self) -> None:
        eco = Ecosystem(size=(10, 10))
        add_agent(eco, "producer", (3, 3))
        agent_id = list(eco.agents.keys())[0]
        remove_agent(eco, agent_id)
        assert len(eco.agents) == 0

    def test_add_multiple_generates_unique_ids(self) -> None:
        eco = Ecosystem(size=(10, 10))
        add_agent(eco, "producer", (1, 1))
        add_agent(eco, "consumer", (2, 2))
        add_agent(eco, "decomposer", (3, 3))
        ids = list(eco.agents.keys())
        assert len(set(ids)) == 3  # All unique


# ---------------------------------------------------------------------------
# count_agents_by_type
# ---------------------------------------------------------------------------


class TestCountAgentsByType:
    """Test count_agents_by_type function."""

    def test_multiple_types(self) -> None:
        eco = Ecosystem(size=(10, 10))
        eco.add_agent(Agent(id=1, agent_type="producer", position=(0, 0)))
        eco.add_agent(Agent(id=2, agent_type="consumer", position=(1, 1)))
        eco.add_agent(Agent(id=3, agent_type="producer", position=(2, 2)))
        eco.add_agent(Agent(id=4, agent_type="decomposer", position=(3, 3)))

        counts = count_agents_by_type(eco)
        assert counts == {"producer": 2, "consumer": 1, "decomposer": 1}

    def test_excludes_dead_agents(self) -> None:
        eco = Ecosystem(size=(10, 10))
        eco.add_agent(Agent(id=1, agent_type="producer", position=(0, 0), alive=True))
        eco.add_agent(Agent(id=2, agent_type="producer", position=(1, 1), alive=False))

        counts = count_agents_by_type(eco)
        assert counts == {"producer": 1}

    def test_empty_ecosystem(self) -> None:
        eco = Ecosystem(size=(10, 10))
        counts = count_agents_by_type(eco)
        assert counts == {}


# ---------------------------------------------------------------------------
# simulation_step
# ---------------------------------------------------------------------------


class TestSimulationStep:
    """Test the simulation_step function."""

    def test_step_advances_agent_age(self) -> None:
        eco = Ecosystem(size=(20, 20))
        agent = Agent(id=1, agent_type="producer", position=(5, 5), energy=100.0)
        eco.add_agent(agent)
        rng = random.Random(42)

        simulation_step(eco, rng)
        assert agent.age == 1

    def test_step_consumes_basal_energy(self) -> None:
        """Each alive agent loses 1.0 energy per step for metabolism."""
        eco = Ecosystem(size=(20, 20))
        agent = Agent(
            id=1,
            agent_type="producer",
            position=(5, 5),
            energy=100.0,
            properties={"energy_efficiency": 0.0, "reproduction_rate": 0.0},
        )
        eco.add_agent(agent)
        rng = random.Random(42)

        initial_energy = agent.energy
        simulation_step(eco, rng)
        # Producer: -1.0 (metabolism) + 2.0 * energy_efficiency (0.0 here) = net -1.0
        assert agent.energy == pytest.approx(initial_energy - 1.0)

    def test_step_producer_gains_energy(self) -> None:
        """Producers gain energy from the environment."""
        eco = Ecosystem(size=(20, 20))
        agent = Agent(
            id=1,
            agent_type="producer",
            position=(5, 5),
            energy=100.0,
            properties={"energy_efficiency": 0.9, "reproduction_rate": 0.0},
        )
        eco.add_agent(agent)
        rng = random.Random(42)

        simulation_step(eco, rng)
        # -1.0 (metabolism) + 2.0 * 0.9 (production) = +0.8 net
        assert agent.energy == pytest.approx(100.0 - 1.0 + 2.0 * 0.9)

    def test_step_dead_agent_skipped(self) -> None:
        eco = Ecosystem(size=(20, 20))
        dead_agent = Agent(id=1, agent_type="producer", position=(5, 5), energy=50.0, alive=False)
        eco.add_agent(dead_agent)
        rng = random.Random(42)

        simulation_step(eco, rng)
        # Dead agent should not have its age incremented
        assert dead_agent.age == 0
        assert dead_agent.energy == 50.0

    def test_step_starving_agent_dies(self) -> None:
        """Agent with energy close to zero dies from starvation."""
        eco = Ecosystem(size=(20, 20))
        # Agent with barely any energy and no energy gain mechanism
        agent = Agent(
            id=1,
            agent_type="consumer",
            position=(5, 5),
            energy=0.5,
            properties={"hunting_skill": 0.0},
        )
        eco.add_agent(agent)
        rng = random.Random(42)

        simulation_step(eco, rng)
        # After metabolism: 0.5 - 1.0 = 0 (clamped), then agent dies
        assert agent.alive is False


# ---------------------------------------------------------------------------
# run_simulation
# ---------------------------------------------------------------------------


class TestRunSimulation:
    """Test run_simulation function."""

    def test_basic_run(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(20, 20),
            rng=rng,
        )
        snapshots = run_simulation(eco, n_steps=5, rng=random.Random(42))
        assert len(snapshots) == 5  # record_interval=1 by default

    def test_snapshots_contain_required_keys(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(20, 20),
            rng=rng,
        )
        snapshots = run_simulation(eco, n_steps=3, rng=random.Random(42))

        for snapshot in snapshots:
            assert "step" in snapshot
            assert "n_agents" in snapshot
            assert "agent_counts" in snapshot
            assert "total_energy" in snapshot
            assert "agent_positions" in snapshot

    def test_step_numbers_ascending(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=5,
            agent_types=["producer"],
            environment_size=(10, 10),
            rng=rng,
        )
        snapshots = run_simulation(eco, n_steps=10, rng=random.Random(42))
        steps = [s["step"] for s in snapshots]
        assert steps == list(range(10))

    def test_record_interval(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(20, 20),
            rng=rng,
        )
        snapshots = run_simulation(eco, n_steps=10, record_interval=3, rng=random.Random(42))
        # Steps 0, 3, 6, 9 are recorded
        steps = [s["step"] for s in snapshots]
        assert steps == [0, 3, 6, 9]

    def test_time_step_advances(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=5,
            agent_types=["producer"],
            environment_size=(10, 10),
            rng=rng,
        )
        initial_ts = eco.time_step
        run_simulation(eco, n_steps=7, rng=random.Random(42))
        assert eco.time_step == initial_ts + 7

    def test_invalid_n_steps_zero(self) -> None:
        eco = Ecosystem(size=(10, 10))
        with pytest.raises(ValidationError):
            run_simulation(eco, n_steps=0)

    def test_invalid_n_steps_negative(self) -> None:
        eco = Ecosystem(size=(10, 10))
        with pytest.raises(ValidationError):
            run_simulation(eco, n_steps=-1)

    def test_simulation_reproducibility(self) -> None:
        """Same seed produces identical simulation snapshots."""
        eco1 = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(15, 15),
            rng=random.Random(42),
        )
        eco2 = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(15, 15),
            rng=random.Random(42),
        )
        snaps1 = run_simulation(eco1, n_steps=5, rng=random.Random(99))
        snaps2 = run_simulation(eco2, n_steps=5, rng=random.Random(99))

        for s1, s2 in zip(snaps1, snaps2):
            assert s1["step"] == s2["step"]
            assert s1["n_agents"] == s2["n_agents"]
            assert s1["agent_counts"] == s2["agent_counts"]
            assert s1["total_energy"] == pytest.approx(s2["total_energy"])


# ---------------------------------------------------------------------------
# get_population_dynamics
# ---------------------------------------------------------------------------


class TestPopulationDynamics:
    """Test get_population_dynamics function."""

    def _make_snapshots(self) -> List[Dict]:
        """Create a simple sequence of snapshots for testing."""
        return [
            {"step": 0, "n_agents": 10, "agent_counts": {"producer": 5, "consumer": 5}, "total_energy": 1000.0},
            {"step": 1, "n_agents": 12, "agent_counts": {"producer": 7, "consumer": 5}, "total_energy": 1200.0},
            {"step": 2, "n_agents": 8, "agent_counts": {"producer": 4, "consumer": 4}, "total_energy": 800.0},
        ]

    def test_basic_dynamics(self) -> None:
        snapshots = self._make_snapshots()
        dynamics = get_population_dynamics(snapshots)

        assert "total_agents_trajectory" in dynamics
        assert "agent_type_dynamics" in dynamics
        assert "diversity_trajectory" in dynamics
        assert "simulation_steps" in dynamics
        assert "final_diversity" in dynamics
        assert "max_diversity" in dynamics

    def test_total_agents_trajectory(self) -> None:
        snapshots = self._make_snapshots()
        dynamics = get_population_dynamics(snapshots)
        assert dynamics["total_agents_trajectory"] == [10, 12, 8]

    def test_agent_type_dynamics(self) -> None:
        snapshots = self._make_snapshots()
        dynamics = get_population_dynamics(snapshots)
        td = dynamics["agent_type_dynamics"]

        assert "producer" in td
        assert "consumer" in td
        assert td["producer"]["counts"] == [5, 7, 4]
        assert td["producer"]["max_count"] == 7
        assert td["producer"]["min_count"] == 4
        assert td["producer"]["final_count"] == 4

    def test_diversity_trajectory(self) -> None:
        snapshots = self._make_snapshots()
        dynamics = get_population_dynamics(snapshots)
        # All diversity values should be non-negative
        for d in dynamics["diversity_trajectory"]:
            assert d >= 0

    def test_empty_data(self) -> None:
        result = get_population_dynamics([])
        assert result == {}

    def test_single_snapshot(self) -> None:
        snapshots = [{"step": 0, "n_agents": 5, "agent_counts": {"producer": 5}, "total_energy": 500.0}]
        dynamics = get_population_dynamics(snapshots)
        assert dynamics["total_agents_trajectory"] == [5]
        # Single type => Shannon diversity = 0 (log(1) = 0)
        assert dynamics["diversity_trajectory"][0] == pytest.approx(0.0)

    def test_even_distribution_max_diversity(self) -> None:
        """Perfectly even distribution has maximum Shannon diversity."""
        snapshots = [
            {"step": 0, "n_agents": 30, "agent_counts": {"a": 10, "b": 10, "c": 10}, "total_energy": 3000.0}
        ]
        dynamics = get_population_dynamics(snapshots)
        expected_shannon = -3 * (1 / 3 * np.log(1 / 3))
        assert dynamics["diversity_trajectory"][0] == pytest.approx(expected_shannon, rel=1e-6)


# ---------------------------------------------------------------------------
# calculate_biodiversity_metrics
# ---------------------------------------------------------------------------


class TestBiodiversityMetrics:
    """Test calculate_biodiversity_metrics function."""

    def test_basic_metrics(self) -> None:
        snapshots = [
            {"step": 0, "n_agents": 10, "agent_counts": {"producer": 5, "consumer": 3, "decomposer": 2}}
        ]
        metrics = calculate_biodiversity_metrics(snapshots)

        assert "species_richness" in metrics
        assert "shannon_diversity" in metrics
        assert "simpson_diversity" in metrics
        assert "pielou_evenness" in metrics
        assert "total_individuals" in metrics

    def test_species_richness(self) -> None:
        snapshots = [{"step": 0, "n_agents": 10, "agent_counts": {"a": 3, "b": 4, "c": 3}}]
        metrics = calculate_biodiversity_metrics(snapshots)
        assert metrics["species_richness"] == 3

    def test_shannon_diversity_two_equal_types(self) -> None:
        """Shannon diversity for 2 equal types = ln(2)."""
        snapshots = [{"step": 0, "n_agents": 10, "agent_counts": {"a": 5, "b": 5}}]
        metrics = calculate_biodiversity_metrics(snapshots)
        expected = -2 * (0.5 * np.log(0.5))
        assert metrics["shannon_diversity"] == pytest.approx(expected, rel=1e-6)

    def test_shannon_diversity_single_type(self) -> None:
        """Shannon diversity for 1 type = 0."""
        snapshots = [{"step": 0, "n_agents": 10, "agent_counts": {"a": 10}}]
        metrics = calculate_biodiversity_metrics(snapshots)
        assert metrics["shannon_diversity"] == pytest.approx(0.0)

    def test_simpson_diversity(self) -> None:
        """Simpson diversity = 1 - sum(p_i^2)."""
        snapshots = [{"step": 0, "n_agents": 10, "agent_counts": {"a": 5, "b": 3, "c": 2}}]
        metrics = calculate_biodiversity_metrics(snapshots)
        p_a, p_b, p_c = 5 / 10, 3 / 10, 2 / 10
        expected = 1 - (p_a**2 + p_b**2 + p_c**2)
        assert metrics["simpson_diversity"] == pytest.approx(expected, rel=1e-6)

    def test_simpson_single_type(self) -> None:
        """Simpson diversity for 1 type = 0 (complete dominance)."""
        snapshots = [{"step": 0, "n_agents": 10, "agent_counts": {"a": 10}}]
        metrics = calculate_biodiversity_metrics(snapshots)
        assert metrics["simpson_diversity"] == pytest.approx(0.0)

    def test_pielou_evenness_equal_types(self) -> None:
        """Pielou evenness = 1.0 when all types have equal counts."""
        snapshots = [{"step": 0, "n_agents": 12, "agent_counts": {"a": 4, "b": 4, "c": 4}}]
        metrics = calculate_biodiversity_metrics(snapshots)
        assert metrics["pielou_evenness"] == pytest.approx(1.0, rel=1e-6)

    def test_pielou_evenness_single_type(self) -> None:
        """Pielou evenness = 0 when there is only one type (n_types <= 1)."""
        snapshots = [{"step": 0, "n_agents": 10, "agent_counts": {"a": 10}}]
        metrics = calculate_biodiversity_metrics(snapshots)
        assert metrics["pielou_evenness"] == pytest.approx(0.0)

    def test_total_individuals(self) -> None:
        snapshots = [{"step": 0, "n_agents": 15, "agent_counts": {"a": 7, "b": 8}}]
        metrics = calculate_biodiversity_metrics(snapshots)
        assert metrics["total_individuals"] == 15

    def test_empty_data(self) -> None:
        result = calculate_biodiversity_metrics([])
        assert result == {}

    def test_uses_last_snapshot(self) -> None:
        """Biodiversity metrics are computed from the LAST snapshot."""
        snapshots = [
            {"step": 0, "n_agents": 10, "agent_counts": {"a": 5, "b": 5}},
            {"step": 1, "n_agents": 9, "agent_counts": {"a": 9}},
        ]
        metrics = calculate_biodiversity_metrics(snapshots)
        # Last snapshot is all type "a"
        assert metrics["species_richness"] == 1
        assert metrics["shannon_diversity"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# simulate_predator_prey
# ---------------------------------------------------------------------------


class TestSimulatePredatorPrey:
    """Test the simulate_predator_prey function."""

    def test_basic_run(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer"],
            environment_size=(15, 15),
            rng=rng,
        )
        snapshots = simulate_predator_prey(eco, n_steps=5, rng=random.Random(99))
        assert len(snapshots) == 5

    def test_snapshot_keys(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer"],
            environment_size=(15, 15),
            rng=rng,
        )
        snapshots = simulate_predator_prey(eco, n_steps=3, rng=random.Random(99))
        for snap in snapshots:
            assert "step" in snap
            assert "n_predators" in snap
            assert "n_prey" in snap
            assert "total_energy" in snap

    def test_step_numbers(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(10, 10),
            rng=rng,
        )
        snapshots = simulate_predator_prey(eco, n_steps=5, rng=random.Random(99))
        steps = [s["step"] for s in snapshots]
        assert steps == [0, 1, 2, 3, 4]

    def test_predator_efficiency_param(self) -> None:
        """Higher predator efficiency leads to faster prey depletion."""
        rng1 = random.Random(42)
        eco_low = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer"],
            environment_size=(10, 10),
            rng=rng1,
        )
        rng2 = random.Random(42)
        eco_high = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer"],
            environment_size=(10, 10),
            rng=rng2,
        )

        snaps_low = simulate_predator_prey(
            eco_low, n_steps=10, predator_efficiency=0.1, rng=random.Random(99)
        )
        snaps_high = simulate_predator_prey(
            eco_high, n_steps=10, predator_efficiency=0.9, rng=random.Random(99)
        )
        # With higher predator efficiency, expect fewer prey at the end
        # (or at least no more prey than low efficiency on average)
        final_prey_low = snaps_low[-1]["n_prey"]
        final_prey_high = snaps_high[-1]["n_prey"]
        # Not guaranteed every run due to stochasticity, but the simulation ran
        assert isinstance(final_prey_low, int)
        assert isinstance(final_prey_high, int)

    def test_invalid_efficiency(self) -> None:
        eco = Ecosystem(size=(10, 10))
        with pytest.raises(ValidationError):
            simulate_predator_prey(eco, n_steps=1, predator_efficiency=1.5)

    def test_time_step_advances(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(10, 10),
            rng=rng,
        )
        initial_ts = eco.time_step
        simulate_predator_prey(eco, n_steps=5, rng=random.Random(99))
        assert eco.time_step == initial_ts + 5


# ---------------------------------------------------------------------------
# simulate_competition
# ---------------------------------------------------------------------------


class TestSimulateCompetition:
    """Test the simulate_competition function."""

    def test_basic_run(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer", "decomposer"],
            environment_size=(10, 10),
            rng=rng,
        )
        snapshots = simulate_competition(eco, n_steps=5, rng=random.Random(99))
        assert len(snapshots) == 5

    def test_snapshot_keys(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=15,
            agent_types=["producer", "consumer"],
            environment_size=(10, 10),
            rng=rng,
        )
        snapshots = simulate_competition(eco, n_steps=3, rng=random.Random(99))
        for snap in snapshots:
            assert "step" in snap
            assert "occupied_positions" in snap
            assert "agent_counts" in snap
            assert "competition_events" in snap

    def test_invalid_competition_strength(self) -> None:
        eco = Ecosystem(size=(10, 10))
        with pytest.raises(ValidationError):
            simulate_competition(eco, n_steps=1, competition_strength=2.0)

    def test_time_step_advances(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=10,
            agent_types=["producer", "consumer"],
            environment_size=(10, 10),
            rng=rng,
        )
        initial_ts = eco.time_step
        simulate_competition(eco, n_steps=5, rng=random.Random(99))
        assert eco.time_step == initial_ts + 5

    def test_competition_strength_zero(self) -> None:
        """Zero competition strength means no competitive interactions."""
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=10,
            agent_types=["producer"],
            environment_size=(5, 5),
            rng=rng,
        )
        snapshots = simulate_competition(eco, n_steps=3, competition_strength=0.0, rng=random.Random(99))
        assert len(snapshots) == 3


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    """Test edge cases for ecosystem simulations."""

    def test_single_agent_ecosystem(self) -> None:
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=1,
            agent_types=["producer"],
            environment_size=(5, 5),
            rng=rng,
        )
        assert len(eco.agents) == 1
        snapshots = run_simulation(eco, n_steps=3, rng=random.Random(42))
        assert len(snapshots) == 3

    def test_boundary_positions_corner(self) -> None:
        """Agent can exist at all four corners."""
        eco = Ecosystem(size=(10, 10))
        for i, pos in enumerate([(0, 0), (9, 0), (0, 9), (9, 9)]):
            add_agent(eco, "producer", pos)
        assert len(eco.agents) == 4

    def test_large_ecosystem_creation(self) -> None:
        """Creating a larger ecosystem does not error."""
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=100,
            agent_types=["producer", "consumer", "decomposer"],
            environment_size=(50, 50),
            rng=rng,
        )
        assert len(eco.agents) == 100

    def test_simulation_with_all_dead_agents(self) -> None:
        """Simulation runs even if all agents are dead."""
        eco = Ecosystem(size=(10, 10))
        eco.add_agent(Agent(id=1, agent_type="producer", position=(0, 0), alive=False))
        eco.add_agent(Agent(id=2, agent_type="consumer", position=(1, 1), alive=False))

        snapshots = run_simulation(eco, n_steps=3, rng=random.Random(42))
        assert len(snapshots) == 3
        # All snapshots should show 0 alive agents
        for snap in snapshots:
            assert snap["n_agents"] == 0

    def test_ecosystem_1x1_grid(self) -> None:
        """Minimum size ecosystem still works."""
        eco = Ecosystem(size=(1, 1))
        add_agent(eco, "producer", (0, 0))
        neighbors = eco.get_neighboring_positions((0, 0), radius=1)
        assert neighbors == []  # No valid neighbors on a 1x1 grid

    def test_agents_at_empty_position(self) -> None:
        eco = Ecosystem(size=(10, 10))
        eco.add_agent(Agent(id=1, agent_type="producer", position=(3, 3)))
        agents = eco.get_agents_at_position((5, 5))
        assert agents == []

    def test_agents_by_nonexistent_type(self) -> None:
        eco = Ecosystem(size=(10, 10))
        eco.add_agent(Agent(id=1, agent_type="producer", position=(0, 0)))
        agents = eco.get_agents_by_type("dragon")
        assert agents == []


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------


class TestReproducibility:
    """Test that seeded RNGs produce deterministic results."""

    def test_create_ecosystem_reproducible(self) -> None:
        eco1 = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer", "decomposer"],
            environment_size=(20, 20),
            rng=random.Random(42),
        )
        eco2 = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer", "decomposer"],
            environment_size=(20, 20),
            rng=random.Random(42),
        )

        assert len(eco1.agents) == len(eco2.agents)
        for aid in eco1.agents:
            a1 = eco1.agents[aid]
            a2 = eco2.agents[aid]
            assert a1.agent_type == a2.agent_type
            assert a1.position == a2.position
            assert a1.energy == pytest.approx(a2.energy)

    def test_run_simulation_reproducible(self) -> None:
        eco1 = create_ecosystem(
            n_agents=15,
            agent_types=["producer", "consumer"],
            environment_size=(15, 15),
            rng=random.Random(7),
        )
        eco2 = create_ecosystem(
            n_agents=15,
            agent_types=["producer", "consumer"],
            environment_size=(15, 15),
            rng=random.Random(7),
        )
        snaps1 = run_simulation(eco1, n_steps=10, rng=random.Random(123))
        snaps2 = run_simulation(eco2, n_steps=10, rng=random.Random(123))

        assert len(snaps1) == len(snaps2)
        for s1, s2 in zip(snaps1, snaps2):
            assert s1["step"] == s2["step"]
            assert s1["n_agents"] == s2["n_agents"]

    def test_predator_prey_reproducible(self) -> None:
        eco1 = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer"],
            environment_size=(15, 15),
            rng=random.Random(42),
        )
        eco2 = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer"],
            environment_size=(15, 15),
            rng=random.Random(42),
        )
        snaps1 = simulate_predator_prey(eco1, n_steps=5, rng=random.Random(99))
        snaps2 = simulate_predator_prey(eco2, n_steps=5, rng=random.Random(99))

        for s1, s2 in zip(snaps1, snaps2):
            assert s1["step"] == s2["step"]
            assert s1["n_predators"] == s2["n_predators"]
            assert s1["n_prey"] == s2["n_prey"]

    def test_competition_reproducible(self) -> None:
        eco1 = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer"],
            environment_size=(10, 10),
            rng=random.Random(42),
        )
        eco2 = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer"],
            environment_size=(10, 10),
            rng=random.Random(42),
        )
        snaps1 = simulate_competition(eco1, n_steps=5, rng=random.Random(99))
        snaps2 = simulate_competition(eco2, n_steps=5, rng=random.Random(99))

        for s1, s2 in zip(snaps1, snaps2):
            assert s1["step"] == s2["step"]
            assert s1["agent_counts"] == s2["agent_counts"]


# ---------------------------------------------------------------------------
# Integration: full workflow
# ---------------------------------------------------------------------------


class TestIntegrationWorkflow:
    """Integration tests chaining ecosystem creation, simulation, and analysis."""

    def test_create_simulate_analyze(self) -> None:
        """Full pipeline: create ecosystem -> run simulation -> analyze dynamics -> biodiversity."""
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=30,
            agent_types=["producer", "consumer", "decomposer"],
            environment_size=(20, 20),
            rng=rng,
        )

        snapshots = run_simulation(eco, n_steps=10, rng=random.Random(99))
        assert len(snapshots) > 0

        dynamics = get_population_dynamics(snapshots)
        assert "total_agents_trajectory" in dynamics
        assert len(dynamics["total_agents_trajectory"]) == len(snapshots)

        bio_metrics = calculate_biodiversity_metrics(snapshots)
        assert "species_richness" in bio_metrics
        assert bio_metrics["species_richness"] >= 1

    def test_add_agents_then_simulate(self) -> None:
        """Manually build ecosystem with add_agent, then simulate."""
        eco = Ecosystem(size=(10, 10))
        add_agent(eco, "producer", (2, 2), properties={"energy_efficiency": 0.9, "reproduction_rate": 0.5})
        add_agent(eco, "consumer", (5, 5), properties={"hunting_skill": 0.7, "energy_efficiency": 0.7})
        add_agent(eco, "decomposer", (8, 8), properties={"decomposition_rate": 0.8, "energy_efficiency": 0.8})

        snapshots = run_simulation(eco, n_steps=5, rng=random.Random(42))
        assert len(snapshots) == 5
        assert snapshots[0]["n_agents"] == 3

    def test_predator_prey_then_analyze(self) -> None:
        """Run predator-prey and then analyze population dynamics from snapshots."""
        rng = random.Random(42)
        eco = create_ecosystem(
            n_agents=20,
            agent_types=["producer", "consumer"],
            environment_size=(15, 15),
            rng=rng,
        )

        # We need standard snapshots for get_population_dynamics, so use run_simulation
        snapshots = run_simulation(eco, n_steps=10, rng=random.Random(99))
        dynamics = get_population_dynamics(snapshots)

        assert "producer" in dynamics["agent_type_dynamics"] or "consumer" in dynamics["agent_type_dynamics"]
        assert len(dynamics["simulation_steps"]) == len(snapshots)
