from dataclasses import dataclass, field
from typing import Dict, Optional

@dataclass(frozen=True)
class Compound:
    """
    Represents a single chemical compound.
    """
    name: str
    formula: Optional[str] = None
    retention_time: Optional[float] = None
    identifiers: Dict[str, str] = field(default_factory=dict)
    
    def __eq__(self, other):
        if not isinstance(other, Compound):
            return False
        # If names match, they are the same (simple heuristic)
        # In rigorous chem, might use InChIKey or similar
        return self.name == other.name and self.retention_time == other.retention_time

    def __hash__(self):
        return hash((self.name, self.retention_time))
