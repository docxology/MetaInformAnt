"""Chemical phenotype analysis module.

Provides tools for analyzing chemical profiles, including cuticular hydrocarbons (CHCs),
metabolomic data, distance matrices, and marker compound identification.
"""

from .compound import Compound
from .profile import ChemicalProfile, distance_matrix, identify_marker_compounds

__all__ = ["Compound", "ChemicalProfile", "distance_matrix", "identify_marker_compounds"]
