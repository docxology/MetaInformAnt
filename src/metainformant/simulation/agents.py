from __future__ import annotations

import random
from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class Agent:
    x: int
    y: int

    def step(self, world: "GridWorld", rng: random.Random) -> None:
        choices = [(1, 0), (-1, 0), (0, 1), (0, -1)]
        dx, dy = rng.choice(choices)
        self.x = (self.x + dx) % world.width
        self.y = (self.y + dy) % world.height


class GridWorld:
    def __init__(self, width: int, height: int, num_agents: int, rng: random.Random | None = None):
        self.width = max(1, width)
        self.height = max(1, height)
        self.rng = rng or random.Random()
        self.agents: List[Agent] = [
            Agent(self.rng.randrange(self.width), self.rng.randrange(self.height))
            for _ in range(max(0, num_agents))
        ]

    def step(self) -> None:
        for agent in self.agents:
            agent.step(self, self.rng)

    def positions(self) -> List[Tuple[int, int]]:
        return [(a.x, a.y) for a in self.agents]


