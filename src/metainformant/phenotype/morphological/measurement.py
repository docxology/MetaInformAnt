import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, Optional

#: Supported length units, each as an exact whole number of micrometres.
_MICROMETRES_PER_UNIT: Dict[str, int] = {"um": 1, "mm": 1_000, "cm": 10_000, "m": 1_000_000}

#: Other spellings of a supported unit: the micro sign (U+00B5) and the Greek letter mu (U+03BC) are both
#: used for the SI symbol of the micrometre.
_UNIT_ALIASES: Dict[str, str] = {"µm": "um", "μm": "um"}


def _canonical_unit(unit: str) -> Optional[str]:
    """The supported unit ``unit`` denotes, or None when it is not a supported length unit."""
    unit = _UNIT_ALIASES.get(unit, unit)
    return unit if unit in _MICROMETRES_PER_UNIT else None


def _scale(value: float, factor: Fraction) -> float:
    """``value * factor``, correctly rounded once (not once per intermediate step).

    The factors between metric units are powers of ten, which binary floating point cannot hold (0.001 is not
    exactly 0.001), so multiplying through float constants leaves errors such as 1 mm -> 1000.0000000000001 um.
    The product is formed exactly instead. Infinities and NaN pass through, and a result beyond the float range
    overflows to an infinity of the same sign.
    """
    if not math.isfinite(value):
        return value * float(factor)
    try:
        return float(Fraction(value) * factor)
    except OverflowError:
        return math.copysign(math.inf, value)


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
        Convert the value, and its uncertainty, to another length unit.

        Supported units: ``um`` (also written ``µm`` or ``μm``), ``mm``, ``cm`` and ``m``; every
        one converts to every other, exactly (each result is the correctly rounded value, not the product of
        inexact float constants). Unit symbols are case-sensitive. Converting to the measurement's own unit
        returns the measurement itself.

        Args:
            target_unit: Unit to convert to.

        Returns:
            A new Measurement in ``target_unit`` (with the same name and description).

        Raises:
            ValueError: If the measurement's unit or ``target_unit`` is not a supported length unit.
        """
        if self.unit == target_unit:
            return self

        source = _canonical_unit(self.unit)
        target = _canonical_unit(target_unit)
        if source is None or target is None:
            supported = ", ".join(_MICROMETRES_PER_UNIT)
            raise ValueError(
                f"Unsupported unit conversion: {self.unit} to {target_unit} (supported units: {supported})"
            )

        factor = Fraction(_MICROMETRES_PER_UNIT[source], _MICROMETRES_PER_UNIT[target])
        new_val = _scale(self.value, factor)
        new_uncert = None if self.uncertainty is None else _scale(self.uncertainty, factor)

        return Measurement(self.name, new_val, target_unit, new_uncert, self.description)
