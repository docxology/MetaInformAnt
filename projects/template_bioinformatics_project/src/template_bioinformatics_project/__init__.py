"""
template_bioinformatics_project — pipeline library for the MetaInformAnt
standalone bioinformatics template.

All data-processing, statistical-analysis, visualisation, and synthetic-data
logic lives here so the numbered scripts in ``scripts/`` stay thin
orchestrators (argparse + path bootstrap + logging + one delegated call).
"""

__all__ = ["analysis", "common", "processing", "synthetic", "visualization"]
