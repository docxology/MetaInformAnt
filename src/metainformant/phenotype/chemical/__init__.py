"""
Chemical phenotype analysis module.

This module provides tools for analyzing chemical profiles, including cuticular hydrocarbons (CHCs)
and other metabolomic data.
"""

from .compound import Compound
from .profile import ChemicalProfile

__all__ = ["Compound", "ChemicalProfile"]
