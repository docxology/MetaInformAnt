"""GWAS workflow orchestration and configuration management.

This is a re-export module that aggregates the public API from:
- workflow_config: GWASWorkflowConfig, config loading/validation, phenotype loading
- workflow_execution: execute_gwas_workflow, run_gwas, run_multi_trait_gwas

All classes and functions are re-exported here for backward compatibility.
"""

from __future__ import annotations

# Re-export configuration and data loading
from metainformant.gwas.workflow.workflow_config import (
    GWASWorkflowConfig,
    _extract_genotype_matrix,
    _load_phenotypes,
    _load_phenotypes_by_id,
    _normalize_config,
    load_gwas_config,
    validate_gwas_config,
)

# Re-export execution functions
from metainformant.gwas.workflow.workflow_execution import (
    execute_gwas_workflow,
    run_gwas,
    run_multi_trait_gwas,
)

__all__ = [
    # Configuration
    "GWASWorkflowConfig",
    "_normalize_config",
    "load_gwas_config",
    "validate_gwas_config",
    "_extract_genotype_matrix",
    "_load_phenotypes",
    "_load_phenotypes_by_id",
    # Execution
    "execute_gwas_workflow",
    "run_gwas",
    "run_multi_trait_gwas",
]
