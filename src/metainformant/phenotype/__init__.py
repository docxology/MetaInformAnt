"""Phenotype module for MetaInformAnt.

Comprehensive phenotypic trait analysis including:
- Morphological measurements and allometric analysis
- Behavioral sequences, ethograms, and diversity metrics
- Chemical profiles (GC-MS, CHC) with distance and marker detection
- Acoustic signals with spectral and temporal analysis
- Electronic tracking (RFID, video, GPS) with movement ecology
- Life course trajectory analysis
- Cross-omic integration (phenotype-genotype, trait-expression)
- Configurable analysis pipelines

Submodules:
    analysis: Life course and temporal analysis.
    behavior: Behavioral phenotype analysis.
    chemical: Chemical phenotype analysis.
    electronic: Electronic/Sensor phenotype analysis.
    morphological: Morphometric phenotype analysis.
    sonic: Acoustic phenotype analysis.
    data: Data loading utilities.
    visualization: Phenotype visualization tools.
    workflow: Pipeline orchestration.
    integration: Cross-omic integration.
"""

from . import analysis
from . import behavior
from . import chemical
from . import data
from . import electronic
from . import gwas_integration
from . import morphological
from . import sonic
from . import visualization
from . import workflow
from . import integration

# Re-export key classes for convenience
from .behavior import Ethogram, BehaviorSequence
from .chemical import Compound, ChemicalProfile
from .electronic import TrackingPoint, Trajectory
from .morphological import Measurement, MorphometricProfile
from .sonic import AcousticSignal
from .workflow import PhenotypePipeline, PipelineConfig, PipelineResult

# GWAS integration imports
from .gwas_integration.phewas import (
    run_phewas,
    phenotype_correlation_matrix,
    genetic_risk_score,
    phenotype_heritability_screen,
    categorize_phenotypes,
)

__all__ = [
    # Submodules
    "analysis",
    "behavior",
    "chemical",
    "data",
    "electronic",
    "morphological",
    "sonic",
    "visualization",
    "workflow",
    "integration",
    # Key classes
    "Ethogram",
    "BehaviorSequence",
    "Compound",
    "ChemicalProfile",
    "TrackingPoint",
    "Trajectory",
    "Measurement",
    "MorphometricProfile",
    "AcousticSignal",
    "PhenotypePipeline",
    "PipelineConfig",
    "PipelineResult",
    # GWAS integration
    "gwas_integration",
    "run_phewas",
    "phenotype_correlation_matrix",
    "genetic_risk_score",
    "phenotype_heritability_screen",
    "categorize_phenotypes",
]
