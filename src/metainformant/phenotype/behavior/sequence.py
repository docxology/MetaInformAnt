from typing import List, Tuple, Dict, Any
from dataclasses import dataclass
from .ethogram import Ethogram
from metainformant.core.utils.errors import ValidationError

@dataclass
class BehaviorEvent:
    """A single occurrence of a behavior."""
    timestamp: float
    code: str
    duration: float = 0.0
    modifiers: Dict[str, Any] = None

    def __post_init__(self):
        if self.modifiers is None:
            self.modifiers = {}

class BehaviorSequence:
    """
    A sequence of behavioral events.
    
    Attributes:
        events (List[BehaviorEvent]): Time-ordered list of events.
        ethogram (Ethogram): The reference ethogram.
    """
    
    def __init__(self, events: List[BehaviorEvent], ethogram: Ethogram):
        self.events = sorted(events, key=lambda x: x.timestamp)
        self.ethogram = ethogram
        self._validate_events()
        
    def _validate_events(self):
        """Ensure all events use codes valid in the ethogram."""
        for event in self.events:
            if not self.ethogram.validate(event.code):
                raise ValidationError(f"Behavior code '{event.code}' not found in ethogram.")
                
    def calculate_time_budget(self) -> Dict[str, float]:
        """
        Calculate the proportion of total time spent in each behavior.
        Assumes events are states with duration.
        """
        total_duration = sum(e.duration for e in self.events)
        if total_duration == 0:
            return {}
            
        budget = {}
        for event in self.events:
            if event.code not in budget:
                budget[event.code] = 0.0
            budget[event.code] += event.duration
            
        return {code: duration / total_duration for code, duration in budget.items()}
    
    def transition_matrix(self) -> Dict[str, Dict[str, float]]:
        """
        Calculate transition probabilities between behaviors.
        """
        transitions = {}
        
        for i in range(len(self.events) - 1):
            current = self.events[i].code
            next_code = self.events[i+1].code
            
            if current not in transitions:
                transitions[current] = {}
            
            if next_code not in transitions[current]:
                transitions[current][next_code] = 0
                
            transitions[current][next_code] += 1
            
        # Normalize
        matrix = {}
        for start, targets in transitions.items():
            total = sum(targets.values())
            matrix[start] = {target: count / total for target, count in targets.items()}
            
        return matrix
