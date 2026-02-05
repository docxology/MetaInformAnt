from dataclasses import dataclass
from typing import Dict, List, Optional

from metainformant.core.utils.errors import ValidationError


@dataclass
class BehaviorDefinition:
    """Definition of a single behavioral act."""

    code: str
    name: str
    description: str
    category: str = "general"


class Ethogram:
    """
    Defines a controlled vocabulary of behaviors (ethogram).

    Attributes:
        behaviors (Dict[str, BehaviorDefinition]): Mapping of behavior codes to definitions.
    """

    def __init__(self, behaviors: Dict[str, str] | Dict[str, BehaviorDefinition]):
        """
        Initialize the ethogram.

        Args:
            behaviors: Dict where key is code. Value can be description string or BehaviorDefinition.
        """
        self.behaviors: Dict[str, BehaviorDefinition] = {}

        for code, value in behaviors.items():
            if isinstance(value, str):
                self.behaviors[code] = BehaviorDefinition(code=code, name=code, description=value)
            elif isinstance(value, BehaviorDefinition):
                self.behaviors[code] = value
            else:
                raise ValidationError(f"Invalid behavior definition for code {code}")

    def validate(self, code: str) -> bool:
        """Check if a behavior code exists in the ethogram."""
        return code in self.behaviors

    def get(self, code: str) -> Optional[BehaviorDefinition]:
        """Get definition for a behavior code."""
        return self.behaviors.get(code)

    def __len__(self) -> int:
        return len(self.behaviors)
