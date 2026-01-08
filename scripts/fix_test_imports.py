#!/usr/bin/env python3
"""Script to fix test imports after module refactoring."""

import os
import re
from pathlib import Path

# Mapping from old import paths to new paths
IMPORT_FIXES = {
    # Visualization
    "metainformant.visualization.genomics.genomics.trees": "metainformant.visualization.genomics.trees",
    # Added mappings for remaining issues
    # Protein structure double-nesting fix
    "metainformant.protein.structure.general.general.general.general.general.alphafold": "metainformant.protein.structure.alphafold",
    "metainformant.protein.structure.general.general.general.general.general.pdb": "metainformant.protein.structure.pdb",
    # RNA amalgkit exports (functions are in amalgkit.py)
    "metainformant.rna.amalgkit": "metainformant.rna.amalgkit.amalgkit",
    "metainformant.rna.amalgkit.build_amalgkit_command": "metainformant.rna.amalgkit.amalgkit.build_amalgkit_command",
    "metainformant.rna.amalgkit.check_cli_available": "metainformant.rna.amalgkit.amalgkit.check_cli_available",
    "metainformant.rna.amalgkit.ensure_cli_available": "metainformant.rna.amalgkit.amalgkit.ensure_cli_available",
    "metainformant.rna.amalgkit.run_amalgkit": "metainformant.rna.amalgkit.amalgkit.run_amalgkit",
    # Math coalescent mapping
    "metainformant.math.coalescent": "metainformant.math.population_genetics.coalescent",

    "metainformant.visualization.genomics.genomics.networks": "metainformant.visualization.genomics.networks",
    "metainformant.visualization.genomics.genomics.expression": "metainformant.visualization.genomics.expression",
    "metainformant.visualization.genomics.genomics.genomics": "metainformant.visualization.genomics.genomics",
    "metainformant.visualization.genomics.genomics": "metainformant.visualization.genomics",
    "metainformant.gwas.visualization.general.comparison": "metainformant.gwas.visualization.comparison",
    "metainformant.gwas.visualization.general.effects": "metainformant.gwas.visualization.effects",
    "metainformant.gwas.visualization.general.genome": "metainformant.gwas.visualization.genome",
    "metainformant.gwas.visualization.general.population": "metainformant.gwas.visualization.population",
    "metainformant.gwas.visualization.general.regional": "metainformant.gwas.visualization.regional",
    "metainformant.gwas.visualization.general.statistical": "metainformant.gwas.visualization.statistical",
    "metainformant.gwas.visualization.general.suite": "metainformant.gwas.visualization.suite",
    "metainformant.gwas.visualization.general.variants": "metainformant.gwas.visualization.variants",
    "metainformant.gwas.visualization.general.general": "metainformant.gwas.visualization.general",
    "metainformant.gwas.visualization.visualization": "metainformant.gwas.visualization.general",
    "metainformant.visualization.basic": "metainformant.visualization.plots.basic",
    "metainformant.visualization.animations": "metainformant.visualization.plots.animations",
    "metainformant.visualization.multidim": "metainformant.visualization.plots.multidim",
    "metainformant.visualization.specialized": "metainformant.visualization.plots.specialized",
    "metainformant.visualization.general": "metainformant.visualization.plots.general",
    "metainformant.visualization.dimred": "metainformant.visualization.analysis.dimred",
    "metainformant.visualization.statistical": "metainformant.visualization.analysis.statistical",
    "metainformant.visualization.quality": "metainformant.visualization.analysis.quality",
    "metainformant.visualization.timeseries": "metainformant.visualization.analysis.timeseries",
    "metainformant.visualization.information": "metainformant.visualization.analysis.information",
    "metainformant.visualization.expression": "metainformant.visualization.genomics.expression",
    "metainformant.visualization.networks": "metainformant.visualization.genomics.networks",
    "metainformant.visualization.trees": "metainformant.visualization.genomics.trees",

    # GWAS
    "metainformant.gwas.visualization": "metainformant.gwas.visualization.general",
    "metainformant.gwas.visualization_comparison": "metainformant.gwas.visualization.comparison",
    "metainformant.gwas.visualization_effects": "metainformant.gwas.visualization.effects",
    "metainformant.gwas.visualization_genome": "metainformant.gwas.visualization.genome",
    "metainformant.gwas.visualization_population": "metainformant.gwas.visualization.population",
    "metainformant.gwas.visualization_regional": "metainformant.gwas.visualization.regional",
    "metainformant.gwas.visualization_statistical": "metainformant.gwas.visualization.statistical",
    "metainformant.gwas.visualization_suite": "metainformant.gwas.visualization.suite",
    "metainformant.gwas.visualization_variants": "metainformant.gwas.visualization.variants",
    "metainformant.gwas.association": "metainformant.gwas.analysis.association",
    "metainformant.gwas.calling": "metainformant.gwas.analysis.calling",
    "metainformant.gwas.correction": "metainformant.gwas.analysis.correction",
    "metainformant.gwas.structure": "metainformant.gwas.analysis.structure",
    "metainformant.gwas.quality": "metainformant.gwas.analysis.quality",
    "metainformant.gwas.download": "metainformant.gwas.data.download",
    "metainformant.gwas.sra_download": "metainformant.gwas.data.sra_download",
    "metainformant.gwas.config": "metainformant.gwas.data.config",
    "metainformant.gwas.workflow": "metainformant.gwas.workflow.workflow",

    # Protein
    "metainformant.protein.alphafold": "metainformant.protein.structure.alphafold",
    "metainformant.protein.pdb": "metainformant.protein.structure.pdb",
    "metainformant.protein.structure_analysis": "metainformant.protein.structure.analysis",
    "metainformant.protein.structure_io": "metainformant.protein.structure.io",
    "metainformant.protein.structure": "metainformant.protein.structure.general",
    "metainformant.protein.secondary": "metainformant.protein.structure.secondary",
    "metainformant.protein.contacts": "metainformant.protein.structure.contacts",
    "metainformant.protein.alignment": "metainformant.protein.sequence.alignment",
    "metainformant.protein.proteomes": "metainformant.protein.sequence.proteomes",
    "metainformant.protein.sequences": "metainformant.protein.sequence.sequences",
    "metainformant.protein.interpro": "metainformant.protein.database.interpro",
    "metainformant.protein.uniprot": "metainformant.protein.database.uniprot",
    "metainformant.protein.visualization": "metainformant.protein.visualization.general",

    # Core
    "metainformant.core.cache": "metainformant.core.io.cache",
    "metainformant.core.disk": "metainformant.core.io.disk",
    "metainformant.core.download": "metainformant.core.io.download",
    "metainformant.core.paths": "metainformant.core.io.paths",
    "metainformant.core.io.io.download": "metainformant.core.io.download",
    "metainformant.core.io.io.io": "metainformant.core.io.io",
    "metainformant.core.io.io.paths": "metainformant.core.io.paths",
    "metainformant.core.io.io.cache": "metainformant.core.io.cache",
    "metainformant.core.io.io.disk": "metainformant.core.io.disk",
    "metainformant.core.config": "metainformant.core.utils.config",
    "metainformant.core.errors": "metainformant.core.utils.errors",
    "metainformant.core.hash": "metainformant.core.utils.hash",
    "metainformant.core.logging": "metainformant.core.utils.logging",
    "metainformant.core.text": "metainformant.core.utils.text",
    "metainformant.core.progress": "metainformant.core.utils.progress",
    "metainformant.core.optional_deps": "metainformant.core.utils.optional_deps",
    "metainformant.core.symbols": "metainformant.core.utils.symbols",
    "metainformant.core.workflow": "metainformant.core.execution.workflow",
    "metainformant.core.parallel": "metainformant.core.execution.parallel",
    "metainformant.core.discovery": "metainformant.core.execution.discovery",
    "metainformant.core.db": "metainformant.core.data.db",
    "metainformant.core.validation": "metainformant.core.data.validation",

    # DNA
    "metainformant.dna.composition": "metainformant.dna.sequence.composition",
    "metainformant.dna.sequence_analysis": "metainformant.dna.sequence.analysis",
    "metainformant.dna.alignment": "metainformant.dna.alignment.pairwise",
    "metainformant.dna.variation": "metainformant.dna.variation.snps",
    "metainformant.dna.expression": "metainformant.dna.expression.degenes",
    "metainformant.dna.population": "metainformant.dna.population.metrics",
    "metainformant.dna.phylogeny": "metainformant.dna.phylogeny.distance",
    "metainformant.dna.io": "metainformant.dna.io.fasta",
    "metainformant.dna.distances": "metainformant.dna.alignment.distances",
    "metainformant.dna.fastq": "metainformant.dna.io.fastq",

    # RNA
    "metainformant.rna.workflow": "metainformant.rna.engine.workflow",
    "metainformant.rna.pipeline": "metainformant.rna.engine.pipeline",
    "metainformant.rna.configs": "metainformant.rna.core.configs",
    "metainformant.rna.deps": "metainformant.rna.core.deps",
    "metainformant.rna.validation": "metainformant.rna.analysis.validation",
    "metainformant.rna.discovery": "metainformant.rna.analysis.discovery",

    # Menu
    "metainformant.menu.display": "metainformant.menu.ui.display",
    "metainformant.menu.navigation": "metainformant.menu.ui.navigation",
    "metainformant.menu.discovery": "metainformant.menu.core.discovery",
    "metainformant.menu.executor": "metainformant.menu.core.executor",

    # Math
    "metainformant.math.coalescent": "metainformant.math.population_genetics.coalescent",
    "metainformant.math.demography": "metainformant.math.population_genetics.demography",
    "metainformant.math.ddm": "metainformant.math.decision_theory.ddm",
    "metainformant.math.epidemiology": "metainformant.math.epidemiology",
    "metainformant.math.evolutionary_dynamics": "metainformant.math.evolutionary_dynamics",

    # Networks
    "metainformant.networks.graph": "metainformant.networks.analysis.graph",
    "metainformant.networks.community": "metainformant.networks.analysis.community",
    "metainformant.networks.pathway": "metainformant.networks.analysis.pathway",
    "metainformant.networks.ppi": "metainformant.networks.interaction.ppi",
    "metainformant.networks.regulatory": "metainformant.networks.interaction.regulatory",

    # Quality
    "metainformant.quality.contamination": "metainformant.quality.analysis.contamination",
    "metainformant.quality.metrics": "metainformant.quality.analysis.metrics",
    "metainformant.quality.fastq": "metainformant.quality.io.fastq",

    # Singlecell
    "metainformant.singlecell.preprocessing": "metainformant.singlecell.data.preprocessing",
    "metainformant.singlecell.integration": "metainformant.singlecell.data.integration",
    "metainformant.singlecell.clustering": "metainformant.singlecell.analysis.clustering",
    "metainformant.singlecell.dimensionality": "metainformant.singlecell.analysis.dimensionality",
    "metainformant.singlecell.trajectory": "metainformant.singlecell.analysis.trajectory",
    "metainformant.singlecell.visualization": "metainformant.singlecell.visualization.visualization",

    # Multiomics
    "metainformant.multiomics.integration": "metainformant.multiomics.analysis.integration",
    "metainformant.multiomics.visualization": "metainformant.multiomics.visualization.visualization",

    # Information
    "metainformant.information.analysis": "metainformant.information.metrics.analysis",
    "metainformant.information.continuous": "metainformant.information.metrics.continuous",
    "metainformant.information.semantic": "metainformant.information.metrics.semantic",
    "metainformant.information.syntactic": "metainformant.information.metrics.syntactic",
    "metainformant.information.estimation": "metainformant.information.metrics.estimation",
    "metainformant.information.integration": "metainformant.information.integration.integration",
    "metainformant.information.networks": "metainformant.information.integration.networks",
    "metainformant.information.workflows": "metainformant.information.workflow.workflows",

    # Life Events
    "metainformant.life_events.models": "metainformant.life_events.models.models",
    "metainformant.life_events.embeddings": "metainformant.life_events.models.embeddings",
    "metainformant.life_events.events": "metainformant.life_events.core.events",
    "metainformant.life_events.config": "metainformant.life_events.core.config",
    "metainformant.life_events.utils": "metainformant.life_events.core.utils",
    "metainformant.life_events.interpretability": "metainformant.life_events.analysis.interpretability",
    "metainformant.life_events.workflow": "metainformant.life_events.workflow.workflow",
    "metainformant.life_events.visualization": "metainformant.life_events.visualization.visualization",

    # Simulation
    "metainformant.simulation.agents": "metainformant.simulation.models.agents",
    "metainformant.simulation.popgen": "metainformant.simulation.models.popgen",
    "metainformant.simulation.rna": "metainformant.simulation.models.rna",
    "metainformant.simulation.sequences": "metainformant.simulation.models.sequences",
    "metainformant.simulation.workflow": "metainformant.simulation.workflow.workflow",
    "metainformant.simulation.visualization": "metainformant.simulation.visualization.visualization",

    # Phenotype
    "metainformant.phenotype.antwiki": "metainformant.phenotype.data.antwiki",
    "metainformant.phenotype.scraper": "metainformant.phenotype.data.scraper",
    "metainformant.phenotype.life_course": "metainformant.phenotype.analysis.life_course",
    "metainformant.phenotype.visualization": "metainformant.phenotype.visualization.visualization",

    # Ecology
    "metainformant.ecology.community": "metainformant.ecology.analysis.community",
    "metainformant.ecology.visualization": "metainformant.ecology.visualization.visualization",

    # Epigenome
    "metainformant.epigenome.atacseq": "metainformant.epigenome.assays.atacseq",
    "metainformant.epigenome.chipseq": "metainformant.epigenome.assays.chipseq",
    "metainformant.epigenome.methylation": "metainformant.epigenome.assays.methylation",
    "metainformant.epigenome.tracks": "metainformant.epigenome.analysis.tracks",
    "metainformant.epigenome.workflow": "metainformant.epigenome.workflow.workflow",
    "metainformant.epigenome.visualization": "metainformant.epigenome.visualization.visualization",

    # ML
    "metainformant.ml.classification": "metainformant.ml.models.classification",
    "metainformant.ml.regression": "metainformant.ml.models.regression",
    "metainformant.ml.features": "metainformant.ml.features.features",
    "metainformant.ml.dimensionality": "metainformant.ml.features.dimensionality",
    "metainformant.ml.validation": "metainformant.ml.evaluation.validation",
}


def fix_imports_in_file(file_path: Path) -> int:
    """Fix imports in a single file.

    Returns:
        Number of replacements made.
    """
    try:
        content = file_path.read_text()
    except UnicodeDecodeError:
        return 0

    original_content = content
    replacements = 0

    for old_path, new_path in IMPORT_FIXES.items():
        # Match "from old_path import ..." or "import old_path"
        pattern = rf'\b{re.escape(old_path)}\b'
        if re.search(pattern, content):
            content = re.sub(pattern, new_path, content)
            replacements += 1

    if replacements > 0:
        file_path.write_text(content)
        print(f"Updated {file_path}: {replacements} replacement(s)")

    return replacements


def main():
    """Fix all test imports."""
    tests_dir = Path("tests")
    if not tests_dir.exists():
        print("tests/ directory not found")
        return

    total_replacements = 0
    files_updated = 0

    for py_file in tests_dir.rglob("*.py"):
        replacements = fix_imports_in_file(py_file)
        if replacements > 0:
            total_replacements += replacements
            files_updated += 1

    print(f"\nTotal: {files_updated} files updated, {total_replacements} replacements")


if __name__ == "__main__":
    main()
