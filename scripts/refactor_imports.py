import os
import re

MAPPING = {
    # Core
    "metainformant.core.cache": "metainformant.core.io.cache",
    "metainformant.core.disk": "metainformant.core.io.disk",
    "metainformant.core.download": "metainformant.core.io.download",
    "metainformant.core.paths": "metainformant.core.io.paths",
    "metainformant.core.io": "metainformant.core.io.io",
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

    # GWAS
    "metainformant.gwas.association": "metainformant.gwas.analysis.association",
    "metainformant.gwas.calling": "metainformant.gwas.analysis.calling",
    "metainformant.gwas.correction": "metainformant.gwas.analysis.correction",
    "metainformant.gwas.quality": "metainformant.gwas.analysis.quality",
    "metainformant.gwas.structure": "metainformant.gwas.analysis.structure",
    "metainformant.gwas.download": "metainformant.gwas.data.download",
    "metainformant.gwas.sra_download": "metainformant.gwas.data.sra_download",
    "metainformant.gwas.config": "metainformant.gwas.data.config",
    "metainformant.gwas.workflow": "metainformant.gwas.workflow.workflow",
    "metainformant.gwas.visualization.visualization": "metainformant.gwas.visualization.general",

    # Protein
    "metainformant.protein.alphafold": "metainformant.protein.structure.alphafold",
    "metainformant.protein.pdb": "metainformant.protein.structure.pdb",
    "metainformant.protein.secondary": "metainformant.protein.structure.secondary",
    "metainformant.protein.contacts": "metainformant.protein.structure.contacts",
    "metainformant.protein.structure_analysis": "metainformant.protein.structure.analysis",
    "metainformant.protein.structure_io": "metainformant.protein.structure.io",
    "metainformant.protein.structure": "metainformant.protein.structure.general",
    "metainformant.protein.alignment": "metainformant.protein.sequence.alignment",
    "metainformant.protein.proteomes": "metainformant.protein.sequence.proteomes",
    "metainformant.protein.sequences": "metainformant.protein.sequence.sequences",
    "metainformant.protein.interpro": "metainformant.protein.database.interpro",
    "metainformant.protein.uniprot": "metainformant.protein.database.uniprot",
    "metainformant.protein.visualization": "metainformant.protein.visualization.general",

    # Epigenome
    "metainformant.epigenome.atacseq": "metainformant.epigenome.assays.atacseq",
    "metainformant.epigenome.chipseq": "metainformant.epigenome.assays.chipseq",
    "metainformant.epigenome.methylation": "metainformant.epigenome.assays.methylation",
    "metainformant.epigenome.tracks": "metainformant.epigenome.analysis.tracks",
    "metainformant.epigenome.workflow": "metainformant.epigenome.workflow.workflow",
    "metainformant.epigenome.visualization": "metainformant.epigenome.visualization.visualization",

    # Networks
    "metainformant.networks.ppi": "metainformant.networks.interaction.ppi",
    "metainformant.networks.regulatory": "metainformant.networks.interaction.regulatory",
    "metainformant.networks.community": "metainformant.networks.analysis.community",
    "metainformant.networks.pathway": "metainformant.networks.analysis.pathway",
    "metainformant.networks.graph": "metainformant.networks.analysis.graph",

    # ML
    "metainformant.ml.classification": "metainformant.ml.models.classification",
    "metainformant.ml.regression": "metainformant.ml.models.regression",
    "metainformant.ml.features": "metainformant.ml.features.features",
    "metainformant.ml.dimensionality": "metainformant.ml.features.dimensionality",
    "metainformant.ml.validation": "metainformant.ml.evaluation.validation",

    # Simulation
    "metainformant.simulation.agents": "metainformant.simulation.models.agents",
    "metainformant.simulation.popgen": "metainformant.simulation.models.popgen",
    "metainformant.simulation.rna": "metainformant.simulation.models.rna",
    "metainformant.simulation.sequences": "metainformant.simulation.models.sequences",
    "metainformant.simulation.workflow": "metainformant.simulation.workflow.workflow",
    "metainformant.simulation.visualization": "metainformant.simulation.visualization.visualization",

    # RNA
    "metainformant.rna.configs": "metainformant.rna.core.configs",
    "metainformant.rna.deps": "metainformant.rna.core.deps",
    "metainformant.rna.environment": "metainformant.rna.core.environment",
    "metainformant.rna.cleanup": "metainformant.rna.core.cleanup",
    "metainformant.rna.orchestration": "metainformant.rna.engine.orchestration",
    "metainformant.rna.pipeline": "metainformant.rna.engine.pipeline",
    "metainformant.rna.workflow": "metainformant.rna.engine.workflow",
    "metainformant.rna.progress_tracker": "metainformant.rna.engine.progress_tracker",
    "metainformant.rna.monitoring": "metainformant.rna.engine.monitoring",
    "metainformant.rna.discovery": "metainformant.rna.engine.discovery",
    "metainformant.rna.amalgkit": "metainformant.rna.amalgkit.amalgkit",
    "metainformant.rna.metadata_filter": "metainformant.rna.amalgkit.metadata_filter",
    "metainformant.rna.genome_prep": "metainformant.rna.amalgkit.genome_prep",
    "metainformant.rna.validation": "metainformant.rna.analysis.validation",
    "metainformant.rna.protein_integration": "metainformant.rna.analysis.protein_integration",

    # Visualization
    "metainformant.visualization.basic": "metainformant.visualization.plots.basic",
    "metainformant.visualization.plots": "metainformant.visualization.plots.general",
    "metainformant.visualization.animations": "metainformant.visualization.plots.animations",
    "metainformant.visualization.specialized": "metainformant.visualization.plots.specialized",
    "metainformant.visualization.multidim": "metainformant.visualization.plots.multidim",
    "metainformant.visualization.genomics": "metainformant.visualization.genomics.genomics",
    "metainformant.visualization.expression": "metainformant.visualization.genomics.expression",
    "metainformant.visualization.networks": "metainformant.visualization.genomics.networks",
    "metainformant.visualization.trees": "metainformant.visualization.genomics.trees",
    "metainformant.visualization.dimred": "metainformant.visualization.analysis.dimred",
    "metainformant.visualization.statistical": "metainformant.visualization.analysis.statistical",
    "metainformant.visualization.information": "metainformant.visualization.analysis.information",
    "metainformant.visualization.quality": "metainformant.visualization.analysis.quality",
    "metainformant.visualization.timeseries": "metainformant.visualization.analysis.timeseries",

    # Singlecell
    "metainformant.singlecell.clustering": "metainformant.singlecell.analysis.clustering",
    "metainformant.singlecell.dimensionality": "metainformant.singlecell.analysis.dimensionality",
    "metainformant.singlecell.trajectory": "metainformant.singlecell.analysis.trajectory",
    "metainformant.singlecell.integration": "metainformant.singlecell.data.integration",
    "metainformant.singlecell.preprocessing": "metainformant.singlecell.data.preprocessing",
    "metainformant.singlecell.visualization": "metainformant.singlecell.visualization.visualization",

    # Phenotype
    "metainformant.phenotype.antwiki": "metainformant.phenotype.data.antwiki",
    "metainformant.phenotype.scraper": "metainformant.phenotype.data.scraper",
    "metainformant.phenotype.life_course": "metainformant.phenotype.analysis.life_course",
    "metainformant.phenotype.visualization": "metainformant.phenotype.visualization.visualization",

    # Ecology
    "metainformant.ecology.community": "metainformant.ecology.analysis.community",
    "metainformant.ecology.visualization": "metainformant.ecology.visualization.visualization",

    # Multiomics
    "metainformant.multiomics.integration": "metainformant.multiomics.analysis.integration",
    "metainformant.multiomics.visualization": "metainformant.multiomics.visualization.visualization",

    # Quality
    "metainformant.quality.contamination": "metainformant.quality.analysis.contamination",
    "metainformant.quality.metrics": "metainformant.quality.analysis.metrics",
    "metainformant.quality.fastq": "metainformant.quality.io.fastq",

    # Menu
    "metainformant.menu.display": "metainformant.menu.ui.display",
    "metainformant.menu.executor": "metainformant.menu.core.executor",
    "metainformant.menu.navigation": "metainformant.menu.ui.navigation",
    "metainformant.menu.discovery": "metainformant.menu.core.discovery",

    # GWAS Workflow specific fixes
    "metainformant.gwas.workflow.quality": "metainformant.gwas.analysis.quality",
    "metainformant.gwas.workflow.structure": "metainformant.gwas.analysis.structure",
    "metainformant.gwas.workflow.association": "metainformant.gwas.analysis.association",
    "metainformant.gwas.workflow.correction": "metainformant.gwas.analysis.correction",
    "metainformant.gwas.workflow.visualization": "metainformant.gwas.visualization.general",

    # Singlecell specific fixes
    "metainformant.singlecell.analysis.preprocessing": "metainformant.singlecell.data.preprocessing",
    "metainformant.singlecell.visualization.preprocessing": "metainformant.singlecell.data.preprocessing",

    # Core IO specific fixes
    "metainformant.core.io.io": "metainformant.core.io.io", # Self-reference fix if needed
    "metainformant.core.utils.errors": "metainformant.core.utils.errors",
    "metainformant.core.data.errors": "metainformant.core.utils.errors",
    "metainformant.core.data.paths": "metainformant.core.io.paths",
    "metainformant.core.execution.io": "metainformant.core.io.io",
    "metainformant.core.execution.paths": "metainformant.core.io.paths",

    # RNA specific fixes
    "metainformant.rna.engine.configs": "metainformant.rna.core.configs",
    "metainformant.rna.amalgkit.amalgkit.metadata_filter": "metainformant.rna.amalgkit.metadata_filter",
    "metainformant.rna.amalgkit.amalgkit": "metainformant.rna.amalgkit.amalgkit",

    # Networks (re-listed for clarity/correction if needed, though already present)
    "metainformant.networks.graph": "metainformant.networks.analysis.graph",
    "metainformant.networks.pathway": "metainformant.networks.analysis.pathway",
    "metainformant.networks.ppi": "metainformant.networks.interaction.ppi",

    # Life Events
    "metainformant.life_events.models": "metainformant.life_events.models.models",
    "metainformant.life_events.embeddings": "metainformant.life_events.models.embeddings",
    "metainformant.life_events.events": "metainformant.life_events.core.events",
    "metainformant.life_events.config": "metainformant.life_events.core.config",
    "metainformant.life_events.utils": "metainformant.life_events.core.utils",
    "metainformant.life_events.interpretability": "metainformant.life_events.analysis.interpretability",
    "metainformant.life_events.workflow": "metainformant.life_events.workflow.workflow",
    "metainformant.life_events.visualization": "metainformant.life_events.visualization.visualization",

    # Information
    "metainformant.information.continuous": "metainformant.information.metrics.continuous",
    "metainformant.information.semantic": "metainformant.information.metrics.semantic",
    "metainformant.information.syntactic": "metainformant.information.metrics.syntactic",
    "metainformant.information.analysis": "metainformant.information.metrics.analysis",
    "metainformant.information.estimation": "metainformant.information.metrics.estimation",
    "metainformant.information.integration": "metainformant.information.integration.integration",
    "metainformant.information.networks": "metainformant.information.integration.networks",
    "metainformant.information.workflows": "metainformant.information.workflow.workflows",
}

def process_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Sort keys by length descending to avoid partial matches
    sorted_keys = sorted(MAPPING.keys(), key=len, reverse=True)
    
    for old_path in sorted_keys:
        new_path = MAPPING[old_path]
        # Replace explicit imports
        content = content.replace(f"from {old_path}", f"from {new_path}")
        content = content.replace(f"import {old_path}", f"import {new_path}")
        
    if content != original_content:
        print(f"Updating imports in {filepath}")
        with open(filepath, 'w') as f:
            f.write(content)

def main():
    root_dir = "src/metainformant"
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith(".py"):
                process_file(os.path.join(root, file))

if __name__ == "__main__":
    main()
