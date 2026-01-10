"""
Behavioral phenotype analysis module.

This module provides tools for defining ethograms, analyzing behavioral sequences,
and calculating time budgets and transition matrices.
"""

from .ethogram import Ethogram
from .sequence import BehaviorSequence

__all__ = ["Ethogram", "BehaviorSequence"]
