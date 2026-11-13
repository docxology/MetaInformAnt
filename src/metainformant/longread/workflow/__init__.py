"""Long-read workflow orchestration module.

Provides end-to-end pipeline orchestration for long-read sequencing analysis
including QC, assembly, methylation, and structural variant calling pipelines.
Ties together the io, quality, analysis, assembly, and visualization subpackages
into cohesive, configurable workflows.

Submodules:
- orchestrator_core: Pipeline execution engine with dependency resolution
- pipeline_stages: Step builders for the QC, assembly, methylation and SV pipelines
- pipelines: Pre-defined pipeline configurations
- reporting: QC and analysis report generation"""

from __future__ import annotations

from . import orchestrator_core, pipeline_stages, pipelines, reporting

__all__ = ["orchestrator_core", "pipeline_stages", "pipelines", "reporting"]
