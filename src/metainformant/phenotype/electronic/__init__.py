"""Electronic phenotype analysis module.

Provides tools for analyzing data from electronic sensors,
including RFID tags, video tracking, and GPS systems.
"""

from .tracking import TrackingPoint, Trajectory, detect_interactions

__all__ = ["TrackingPoint", "Trajectory", "detect_interactions"]
