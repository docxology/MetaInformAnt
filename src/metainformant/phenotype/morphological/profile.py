from typing import Dict, Optional, Any
from .measurement import Measurement
from metainformant.core.utils.errors import ValidationError


class MorphometricProfile:
    """
    A collection of morphometric measurements for a single specimen.
    """

    def __init__(self, specimen_id: str, measurements: List[Measurement] = None, metadata: Dict[str, Any] = None):
        self.specimen_id = specimen_id
        self._measurements: Dict[str, Measurement] = {}
        if measurements:
            for m in measurements:
                self.add(m)
        self.metadata = metadata or {}

    def add(self, measurement: Measurement):
        """Add or replace a measurement."""
        self._measurements[measurement.name] = measurement

    def get(self, name: str) -> Optional[Measurement]:
        """Retrieve a measurement by name."""
        return self._measurements.get(name)

    def calculate_index(self, name: str, numerator: str, denominator: str) -> float:
        """
        Calculate a morphometric index (ratio * 100).
        Example: Cephalic Index (CI) = (Head Width / Head Length) * 100
        """
        m_num = self.get(numerator)
        m_den = self.get(denominator)

        if not m_num or not m_den:
            raise ValidationError(f"Missing measurements for index {name}: {numerator}, {denominator}")

        # Ensure units match or standarize (simple check for now)
        if m_num.unit != m_den.unit:
            # Basic conversion attempt
            try:
                m_num = m_num.convert(m_den.unit)
            except ValueError:
                raise ValidationError(f"Unit mismatch: {m_num.unit} vs {m_den.unit}")

        if m_den.value == 0:
            raise ValueError("Denominator is zero")

        return (m_num.value / m_den.value) * 100.0
