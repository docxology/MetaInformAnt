from __future__ import annotations

import random
from dataclasses import dataclass
from typing import List, Tuple

from ..core import logging, validation

logger = logging.get_logger(__name__)


@dataclass
class Agent:
    """Agent in a grid world simulation.
    
    Attributes:
        x: X coordinate in grid
        y: Y coordinate in grid
    """
    x: int
    y: int

    def step(self, world: "GridWorld", rng: random.Random) -> None:
        """Move agent randomly one step in grid (with wraparound).
        
        Args:
            world: Grid world environment
            rng: Random number generator
        """
        choices = [(1, 0), (-1, 0), (0, 1), (0, -1)]
        dx, dy = rng.choice(choices)
        self.x = (self.x + dx) % world.width
        self.y = (self.y + dy) % world.height


class GridWorld:
    """Grid-based world for agent-based simulation.
    
    Attributes:
        width: Grid width
        height: Grid height
        rng: Random number generator
        agents: List of agents in the world
    """
    def __init__(self, width: int, height: int, num_agents: int, rng: random.Random | None = None):
        """Initialize grid world with agents.
        
        Args:
            width: Grid width (minimum 1)
            height: Grid height (minimum 1)
            num_agents: Number of agents to create
            rng: Random number generator (default: new Random())
            
        Raises:
            ValidationError: If parameters are invalid
        """
        validation.validate_type(width, int, "width")
        validation.validate_type(height, int, "height")
        validation.validate_type(num_agents, int, "num_agents")
        validation.validate_range(width, min_val=1, name="width")
        validation.validate_range(height, min_val=1, name="height")
        validation.validate_range(num_agents, min_val=0, name="num_agents")
        
        self.width = max(1, width)
        self.height = max(1, height)
        self.rng = rng or random.Random()
        self.agents: List[Agent] = [
            Agent(self.rng.randrange(self.width), self.rng.randrange(self.height)) for _ in range(max(0, num_agents))
        ]
        
        logger.info(f"Initialized GridWorld: {width}x{height} with {num_agents} agents")

    def step(self) -> None:
        """Advance simulation one time step (all agents move)."""
        for agent in self.agents:
            agent.step(self, self.rng)

    def positions(self) -> List[Tuple[int, int]]:
        """Get current positions of all agents.
        
        Returns:
            List of (x, y) tuples for each agent
        """
        return [(a.x, a.y) for a in self.agents]
