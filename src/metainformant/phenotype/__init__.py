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

from . import (
    analysis,
    behavior,
    chemical,
    data,
    electronic,
    gwas_integration,
    integration,
    morphological,
    sonic,
    visualization,
    workflow,
)

# Re-export key classes for convenience
from .behavior import BehaviorSequence, Ethogram
from .chemical import ChemicalProfile, Compound
from .electronic import TrackingPoint, Trajectory

# GWAS integration imports
from .gwas_integration.phewas import (
    categorize_phenotypes,
    genetic_risk_score,
    phenotype_correlation_matrix,
    phenotype_heritability_screen,
    run_phewas,
)
from .morphological import Measurement, MorphometricProfile
from .sonic import AcousticSignal
from .workflow import PhenotypePipeline, PipelineConfig, PipelineResult

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
