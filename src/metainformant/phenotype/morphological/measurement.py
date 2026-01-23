from dataclasses import dataclass
from typing import Optional


@dataclass
class Measurement:
    """
    A single morphometric measurement.
    """

    name: str
    value: float
    unit: str = "mm"
    uncertainty: Optional[float] = None
    description: Optional[str] = None

    def convert(self, target_unit: str) -> "Measurement":
        """
        Simple conversion (placeholder for comprehensive unit lib).
        Currently supports mm <-> cm <-> m.
        """
        if self.unit == target_unit:
            return self

        factors = {"mm": 0.001, "cm": 0.01, "m": 1.0, "um": 0.000001}

        if self.unit not in factors or target_unit not in factors:
            raise ValueError(f"Unsupported unit conversion: {self.unit} to {target_unit}")

        val_meters = self.value * factors[self.unit]
        new_val = val_meters / factors[target_unit]

        new_uncert = None
        if self.uncertainty is not None:
            uncert_meters = self.uncertainty * factors[self.unit]
            new_uncert = uncert_meters / factors[target_unit]

        return Measurement(self.name, new_val, target_unit, new_uncert, self.description)
