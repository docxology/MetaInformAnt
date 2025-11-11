"""Ecology domain functionality."""

from .community import (
    bray_curtis_dissimilarity,
    chao1_estimator,
    community_metrics,
    functional_diversity_metrics,
    jaccard_similarity,
    pielou_evenness,
    rank_abundance_distribution,
    rarefaction_curve,
    shannon_diversity,
    simpson_diversity,
    sorensen_similarity,
    species_accumulation_curve,
    species_richness,
    beta_diversity_partitioning,
)
from .environmental import (
    analyze_environmental_gradient,
    spatial_autocorrelation,
)
from .interactions import (
    build_interaction_network,
    calculate_interaction_strength,
    identify_keystone_species,
)

try:
    from .workflow import run_community_analysis_workflow
    _workflow_available = True
except ImportError:
    _workflow_available = False

__all__ = [
    # Diversity metrics
    "shannon_diversity",
    "simpson_diversity",
    "species_richness",
    "pielou_evenness",
    "chao1_estimator",
    "community_metrics",
    # Beta diversity
    "bray_curtis_dissimilarity",
    "jaccard_similarity",
    "sorensen_similarity",
    "beta_diversity_partitioning",
    # Curves and distributions
    "rarefaction_curve",
    "species_accumulation_curve",
    "rank_abundance_distribution",
    # Functional diversity
    "functional_diversity_metrics",
    # Interactions
    "build_interaction_network",
    "calculate_interaction_strength",
    "identify_keystone_species",
    # Environmental analysis
    "analyze_environmental_gradient",
    "spatial_autocorrelation",
]

if _workflow_available:
    __all__.append("run_community_analysis_workflow")
