"""Shim module for metainformant.math.coalescent.

Re-exports the coalescent functions from the population_genetics subpackage to maintain backward compatibility.
"""

from .population_genetics.coalescent import *  # noqa: F403,F401
