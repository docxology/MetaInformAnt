"""
Electronic phenotype analysis module.

This module provides tools for analyzing data from electronic sensors,
including RFID tags and tracking systems.
"""

from .tracking import TrackingPoint, Trajectory

__all__ = ["TrackingPoint", "Trajectory"]
