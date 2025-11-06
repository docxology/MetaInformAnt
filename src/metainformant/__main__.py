from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from .rna.configs import AmalgkitRunLayout, SpeciesProfile, build_step_params

# defer DNA imports until the DNA subcommand is used to avoid optional deps at import time
from .rna.workflow import AmalgkitWorkflowConfig, execute_workflow, plan_workflow_with_params

# Explicitly import the test runner module to avoid confusion with the repo's root-level tests/
from .tests.runner import run_all_tests


def main() -> None:
    parser = argparse.ArgumentParser(prog="metainformant", description="METAINFORMANT CLI")
    subparsers = parser.add_subparsers(dest="command")

    # setup subcommand
    setup_parser = subparsers.add_parser("setup", help="Run repository setup (uv, deps)")
    setup_parser.add_argument("--with-amalgkit", action="store_true", help="Install AMALGKIT")
    setup_parser.add_argument("--ncbi-email", default="", help="Export NCBI_EMAIL during setup")

    # dna subcommand
    dna_parser = subparsers.add_parser("dna", help="DNA-related operations")
    dna_sub = dna_parser.add_subparsers(dest="dna_cmd")
    fetch = dna_sub.add_parser("fetch", help="Fetch genome by assembly accession")
    fetch.add_argument("--assembly", required=True, help="NCBI assembly accession, e.g., GCF_*")
    
    dna_align = dna_sub.add_parser("align", help="Sequence alignment operations")
    dna_align.add_argument("--input", type=Path, required=True, help="Input FASTA file with sequences")
    dna_align.add_argument("--output", type=Path, default=Path("output/dna/alignment"), help="Output directory")
    dna_align.add_argument("--method", choices=["pairwise", "msa"], default="pairwise", help="Alignment method")
    dna_align.add_argument("--pairwise-type", choices=["global", "local"], default="global", help="Pairwise alignment type")
    
    dna_phylogeny = dna_sub.add_parser("phylogeny", help="Phylogenetic tree construction")
    dna_phylogeny.add_argument("--input", type=Path, required=True, help="Input alignment file or FASTA")
    dna_phylogeny.add_argument("--output", type=Path, default=Path("output/dna/phylogeny"), help="Output directory")
    dna_phylogeny.add_argument("--method", choices=["neighbor_joining"], default="neighbor_joining", help="Tree construction method")
    dna_phylogeny.add_argument("--format", choices=["newick", "json"], default="newick", help="Output format")
    
    dna_population = dna_sub.add_parser("population", help="Population genetics analysis")
    dna_population.add_argument("--input", type=Path, required=True, help="Input sequences file (FASTA)")
    dna_population.add_argument("--output", type=Path, default=Path("output/dna/population"), help="Output directory")
    dna_population.add_argument("--statistics", action="append", choices=["pi", "tajima_d", "fst", "all"], default=["all"], help="Statistics to calculate")
    
    dna_variants = dna_sub.add_parser("variants", help="Variant calling and analysis")
    dna_variants.add_argument("--input", type=Path, required=True, help="Input BAM/VCF file")
    dna_variants.add_argument("--reference", type=Path, help="Reference genome FASTA")
    dna_variants.add_argument("--output", type=Path, default=Path("output/dna/variants"), help="Output directory")
    dna_variants.add_argument("--format", choices=["vcf", "bam"], default="vcf", help="Input format")

    # rna subcommand
    rna_parser = subparsers.add_parser("rna", help="RNA-related operations (amalgkit)")
    rna_sub = rna_parser.add_subparsers(dest="rna_cmd")
    rna_plan = rna_sub.add_parser("plan", help="Plan the amalgkit workflow")
    rna_plan.add_argument("--work-dir", required=True, help="Working directory for the run")
    rna_plan.add_argument("--threads", type=int, default=4)
    rna_plan.add_argument("--species", action="append", default=[], help="Species name (repeatable)")

    rna_run = rna_sub.add_parser("run", help="Execute the amalgkit workflow")
    rna_run_cfg = rna_sub.add_parser("run-config", help="Run amalgkit workflow from a config file")
    rna_run_cfg.add_argument("--config", required=True, help="Path to YAML/TOML/JSON config file under config/")
    rna_run_cfg.add_argument("--check", action="store_true", help="Stop on first failure")
    rna_run.add_argument("--work-dir", required=True, help="Working directory for the run")
    rna_run.add_argument("--threads", type=int, default=4)
    rna_run.add_argument("--species", action="append", default=[], help="Species name (repeatable)")
    rna_run.add_argument("--check", action="store_true", help="Stop on first failure")

    rna_plan_species = rna_sub.add_parser("plan-species", help="Plan workflow with species/tissue params")
    rna_plan_species.add_argument("--work-dir", required=True)
    rna_plan_species.add_argument("--threads", type=int, default=4)
    rna_plan_species.add_argument("--taxon-id", type=int, required=False)
    rna_plan_species.add_argument("--tissue", action="append", default=[])

    # NEW: plan-config subcommand (plan from config without executing)
    rna_plan_cfg = rna_sub.add_parser("plan-config", help="Plan the amalgkit workflow from a config file")
    rna_plan_cfg.add_argument("--config", required=True, help="Path to YAML/TOML/JSON config file under config/")

    # tests subcommand
    tests_parser = subparsers.add_parser("tests", help="Run repository test suite")
    tests_parser.add_argument("pytest_args", nargs=argparse.REMAINDER, help="Arguments passed to pytest")

    # protein subcommand
    protein_parser = subparsers.add_parser("protein", help="Protein-related operations")
    protein_sub = protein_parser.add_subparsers(dest="protein_cmd")
    taxon_ids = protein_sub.add_parser("taxon-ids", help="Print cleaned taxon IDs from a file")
    taxon_ids.add_argument("--file", required=True, help="Path to taxon id list file")
    comp = protein_sub.add_parser("comp", help="Compute amino acid composition for sequences in FASTA")
    comp.add_argument("--fasta", required=True, help="Path to protein FASTA file")
    rmsd_ca = protein_sub.add_parser("rmsd-ca", help="Compute Kabsch RMSD using CA atoms from two PDB files")
    rmsd_ca.add_argument("--pdb-a", required=True)
    rmsd_ca.add_argument("--pdb-b", required=True)
    
    protein_align = protein_sub.add_parser("align", help="Protein sequence alignment")
    protein_align.add_argument("--input", type=Path, required=True, help="Input FASTA file with sequences")
    protein_align.add_argument("--output", type=Path, default=Path("output/protein/alignment"), help="Output directory")
    protein_align.add_argument("--method", choices=["needleman_wunsch", "pairwise"], default="needleman_wunsch", help="Alignment method")
    
    protein_structure = protein_sub.add_parser("structure", help="Protein structure analysis")
    protein_structure.add_argument("--input", type=Path, required=True, help="Input FASTA file or PDB file")
    protein_structure.add_argument("--output", type=Path, default=Path("output/protein/structure"), help="Output directory")
    protein_structure.add_argument("--analyze", choices=["secondary", "stability", "domains", "ptm", "all"], default=["all"], action="append", help="Analysis type")
    
    protein_domains = protein_sub.add_parser("domains", help="Protein domain annotation")
    protein_domains.add_argument("--input", type=Path, required=True, help="Input FASTA file")
    protein_domains.add_argument("--output", type=Path, default=Path("output/protein/domains"), help="Output directory")
    protein_domains.add_argument("--uniprot", action="store_true", help="Fetch InterPro domains from UniProt")
    
    protein_families = protein_sub.add_parser("families", help="Protein family classification")
    protein_families.add_argument("--input", type=Path, required=True, help="Input FASTA file")
    protein_families.add_argument("--output", type=Path, default=Path("output/protein/families"), help="Output directory")
    protein_families.add_argument("--method", choices=["similarity", "pattern"], default="similarity", help="Classification method")

    # math subcommand
    math_parser = subparsers.add_parser("math", help="Math experiments and demos")
    math_sub = math_parser.add_subparsers(dest="math_cmd")
    # selection sub-tree (wire lazily to keep import cost low)
    from .math.selection_experiments.cli import add_math_selection_subparser as _add_sel

    _add_sel(math_sub)
    
    math_popgen = math_sub.add_parser("popgen", help="Population genetics models")
    math_popgen.add_argument("--input", type=Path, help="Input sequence file (FASTA)")
    math_popgen.add_argument("--output", type=Path, default=Path("output/math/popgen"), help="Output directory")
    math_popgen.add_argument("--analysis", choices=["tajima_d", "fst", "ld", "coalescent", "all"], default=["all"], action="append", help="Analysis type")
    
    math_coalescent = math_sub.add_parser("coalescent", help="Coalescent simulations")
    math_coalescent.add_argument("--n-samples", type=int, default=10, help="Number of samples")
    math_coalescent.add_argument("--n-loci", type=int, default=1000, help="Number of loci")
    math_coalescent.add_argument("--output", type=Path, default=Path("output/math/coalescent"), help="Output directory")
    math_coalescent.add_argument("--model", choices=["constant", "exponential", "bottleneck"], default="constant", help="Demographic model")

    # gwas subcommand
    gwas_parser = subparsers.add_parser("gwas", help="GWAS (Genome-Wide Association Studies) workflow")
    gwas_sub = gwas_parser.add_subparsers(dest="gwas_cmd")
    gwas_run = gwas_sub.add_parser("run", help="Run GWAS workflow from configuration file")
    gwas_run.add_argument("--config", required=True, help="Path to GWAS YAML/TOML/JSON config file")
    gwas_run.add_argument("--check", action="store_true", help="Validate configuration only, do not execute")

    # ontology subcommand
    ontology_parser = subparsers.add_parser("ontology", help="Ontology analysis operations (Gene Ontology, functional annotation)")
    ontology_sub = ontology_parser.add_subparsers(dest="ontology_cmd")
    ontology_run = ontology_sub.add_parser(
        "run",
        help="Run ontology analysis workflow",
        description="Execute ontology analysis including GO term queries, enrichment analysis, and ontology summaries",
    )
    ontology_run.add_argument("--go", type=Path, help="Gene Ontology OBO file (required for most operations)")
    ontology_run.add_argument("--output", type=Path, default=Path("output/ontology"), help="Output directory (default: output/ontology)")
    ontology_run.add_argument("--query-term", type=str, help="Query specific GO term ID (e.g., GO:0008150)")
    ontology_run.add_argument("--ancestors", action="store_true", help="Get all ancestor terms (broader terms) for query term")
    ontology_run.add_argument("--descendants", action="store_true", help="Get all descendant terms (more specific terms) for query term")

    # phenotype subcommand
    phenotype_parser = subparsers.add_parser("phenotype", help="Phenotype analysis operations (trait analysis, curation)")
    phenotype_sub = phenotype_parser.add_subparsers(dest="phenotype_cmd")
    phenotype_run = phenotype_sub.add_parser(
        "run",
        help="Run phenotype analysis workflow",
        description="Analyze phenotypic traits including statistics, correlations, and AntWiki data integration",
    )
    phenotype_run.add_argument("--input", type=Path, required=True, help="Input phenotype data file (JSON, CSV, or TSV format)")
    phenotype_run.add_argument("--output", type=Path, default=Path("output/phenotype"), help="Output directory (default: output/phenotype)")
    phenotype_run.add_argument("--analyze-statistics", action="store_true", help="Calculate basic statistics for numeric traits")
    phenotype_run.add_argument("--analyze-correlations", action="store_true", help="Calculate correlations between traits")

    # networks subcommand
    networks_parser = subparsers.add_parser("networks", help="Network analysis operations (PPI, regulatory networks, pathways)")
    networks_sub = networks_parser.add_subparsers(dest="networks_cmd")
    networks_run = networks_sub.add_parser(
        "run",
        help="Run network analysis workflow",
        description="Analyze biological networks including metrics, community detection, and centrality measures",
    )
    networks_run.add_argument("--input", type=Path, required=True, help="Input interaction file (edge list or network format)")
    networks_run.add_argument("--output", type=Path, default=Path("output/networks"), help="Output directory (default: output/networks)")
    networks_run.add_argument("--analyze-metrics", action="store_true", help="Calculate network metrics (density, clustering, etc.)")
    networks_run.add_argument("--detect-communities", action="store_true", help="Detect network communities (Louvain, Leiden algorithms)")
    networks_run.add_argument("--analyze-centrality", action="store_true", help="Calculate centrality measures (degree, betweenness, etc.)")

    # multiomics subcommand
    multiomics_parser = subparsers.add_parser("multiomics", help="Multi-omics integration operations (cross-platform data harmonization)")
    multiomics_sub = multiomics_parser.add_subparsers(dest="multiomics_cmd")
    multiomics_run = multiomics_sub.add_parser(
        "run",
        help="Run multi-omics integration workflow",
        description="Integrate and analyze multiple omics datasets (genomics, transcriptomics, proteomics) using joint dimensionality reduction",
    )
    multiomics_run.add_argument("--genomics", type=Path, help="Genomics data file (variant or genotype data)")
    multiomics_run.add_argument("--transcriptomics", type=Path, help="Transcriptomics data file (expression matrix)")
    multiomics_run.add_argument("--proteomics", type=Path, help="Proteomics data file (protein abundance matrix)")
    multiomics_run.add_argument("--output", type=Path, default=Path("output/multiomics"), help="Output directory (default: output/multiomics)")
    multiomics_run.add_argument("--joint-pca", action="store_true", help="Perform joint Principal Component Analysis across omics")
    multiomics_run.add_argument("--joint-nmf", action="store_true", help="Perform joint Non-negative Matrix Factorization")
    multiomics_run.add_argument("--canonical-correlation", action="store_true", help="Perform Canonical Correlation Analysis (CCA)")

    # singlecell subcommand
    singlecell_parser = subparsers.add_parser("singlecell", help="Single-cell analysis operations (scRNA-seq pipeline)")
    singlecell_sub = singlecell_parser.add_subparsers(dest="singlecell_cmd")
    singlecell_run = singlecell_sub.add_parser(
        "run",
        help="Run single-cell analysis workflow",
        description="Complete scRNA-seq analysis pipeline including QC, normalization, dimensionality reduction, and clustering",
    )
    singlecell_run.add_argument("--input", type=Path, required=True, help="Input count matrix file (genes x cells)")
    singlecell_run.add_argument("--output", type=Path, default=Path("output/singlecell"), help="Output directory (default: output/singlecell)")
    singlecell_run.add_argument("--qc", action="store_true", help="Perform quality control (filter cells/genes, calculate QC metrics)")
    singlecell_run.add_argument("--normalize", action="store_true", help="Normalize counts (log normalization, scaling)")
    singlecell_run.add_argument("--cluster", action="store_true", help="Perform clustering (Leiden, Louvain algorithms) and cell type identification")

    # quality subcommand
    quality_parser = subparsers.add_parser("quality", help="Quality control operations (FASTQ analysis, contamination detection)")
    quality_sub = quality_parser.add_subparsers(dest="quality_cmd")
    quality_run = quality_sub.add_parser(
        "run",
        help="Run quality control workflow",
        description="Comprehensive quality assessment for sequencing data including FASTQ metrics and contamination detection",
    )
    quality_run.add_argument("--fastq", type=Path, help="Input FASTQ file for quality analysis")
    quality_run.add_argument("--output", type=Path, default=Path("output/quality"), help="Output directory (default: output/quality)")
    quality_run.add_argument("--analyze-fastq", action="store_true", help="Analyze FASTQ quality metrics (per-base quality, GC content, duplication)")
    quality_run.add_argument("--detect-contamination", action="store_true", help="Detect contamination (adapter, vector, cross-species, rRNA)")

    # simulation subcommand
    simulation_parser = subparsers.add_parser("simulation", help="Simulation operations (synthetic data generation, agent-based models)")
    simulation_sub = simulation_parser.add_subparsers(dest="simulation_cmd")
    simulation_run = simulation_sub.add_parser(
        "run",
        help="Run simulation workflow",
        description="Generate synthetic biological data and run agent-based simulations for method validation and benchmarking",
    )
    simulation_run.add_argument("--model", type=str, required=True, choices=["sequences", "agents", "expression"], help="Simulation model type")
    simulation_run.add_argument("--output", type=Path, default=Path("output/simulation"), help="Output directory (default: output/simulation)")
    simulation_run.add_argument("--n", type=int, default=100, help="Number of sequences to generate (for sequences model)")
    simulation_run.add_argument("--steps", type=int, default=100, help="Number of simulation steps (for agents model)")

    # visualization subcommand
    visualization_parser = subparsers.add_parser("visualization", help="Visualization operations (plots, heatmaps, animations, trees)")
    visualization_sub = visualization_parser.add_subparsers(dest="visualization_cmd")
    visualization_run = visualization_sub.add_parser(
        "run",
        help="Run visualization workflow",
        description="Generate publication-quality plots and visualizations from biological data",
    )
    visualization_run.add_argument("--input", type=Path, required=True, help="Input data file (CSV, TSV, JSON, or data format)")
    visualization_run.add_argument("--plot-type", type=str, required=True, choices=["lineplot", "heatmap", "animation", "histogram"], help="Type of plot to generate")
    visualization_run.add_argument("--output", type=Path, default=Path("output/visualization"), help="Output directory (default: output/visualization)")

    # epigenome subcommand
    epigenome_parser = subparsers.add_parser("epigenome", help="Epigenome analysis operations (DNA methylation, chromatin tracks)")
    epigenome_sub = epigenome_parser.add_subparsers(dest="epigenome_cmd")
    epigenome_run = epigenome_sub.add_parser(
        "run",
        help="Run epigenome analysis workflow",
        description="Analyze epigenetic modifications including DNA methylation patterns and chromatin accessibility tracks",
    )
    epigenome_run.add_argument("--methylation", type=Path, help="Methylation data file (CpG table format)")
    epigenome_run.add_argument("--bedgraph", type=Path, help="BedGraph track file (chromatin accessibility or ChIP-seq data)")
    epigenome_run.add_argument("--output", type=Path, default=Path("output/epigenome"), help="Output directory (default: output/epigenome)")
    epigenome_run.add_argument("--compute-beta", action="store_true", help="Compute beta values (methylation levels) from CpG counts")

    # ecology subcommand
    ecology_parser = subparsers.add_parser("ecology", help="Ecology analysis operations (community diversity, species abundance)")
    ecology_sub = ecology_parser.add_subparsers(dest="ecology_cmd")
    ecology_run = ecology_sub.add_parser(
        "run",
        help="Run ecology analysis workflow",
        description="Analyze ecological communities including diversity metrics, species richness, and beta diversity",
    )
    ecology_run.add_argument("--input", type=Path, required=True, help="Input abundance table (species x samples matrix)")
    ecology_run.add_argument("--output", type=Path, default=Path("output/ecology"), help="Output directory (default: output/ecology)")
    ecology_run.add_argument("--diversity", action="store_true", help="Calculate diversity indices (Shannon, Simpson, Chao1)")
    ecology_run.add_argument("--beta-diversity", action="store_true", help="Calculate beta diversity (Bray-Curtis, Jaccard dissimilarity)")

    # ml subcommand
    ml_parser = subparsers.add_parser("ml", help="Machine learning operations (classification, regression, feature selection)")
    ml_sub = ml_parser.add_subparsers(dest="ml_cmd")
    ml_run = ml_sub.add_parser(
        "run",
        help="Run ML pipeline workflow",
        description="Complete machine learning pipeline for biological data including feature selection, model training, and validation",
    )
    ml_run.add_argument("--features", type=Path, required=True, help="Feature matrix file (samples x features)")
    ml_run.add_argument("--labels", type=Path, help="Labels file (class labels or continuous outcomes)")
    ml_run.add_argument("--output", type=Path, default=Path("output/ml"), help="Output directory (default: output/ml)")
    ml_run.add_argument("--classify", action="store_true", help="Perform classification (supervised learning for discrete outcomes)")
    ml_run.add_argument("--regress", action="store_true", help="Perform regression (supervised learning for continuous outcomes)")
    ml_run.add_argument("--feature-selection", action="store_true", help="Perform feature selection (identify important features)")

    # information subcommand
    information_parser = subparsers.add_parser("information", help="Information theory operations")
    information_sub = information_parser.add_subparsers(dest="information_cmd")
    information_entropy = information_sub.add_parser("entropy", help="Calculate Shannon entropy")
    information_entropy.add_argument("--input", type=Path, required=True, help="Input sequence file or data file")
    information_entropy.add_argument("--output", type=Path, default=Path("output/information"), help="Output directory")
    information_entropy.add_argument("--k", type=int, default=1, help="K-mer size for sequence entropy")
    information_mi = information_sub.add_parser("mutual-information", help="Calculate mutual information")
    information_mi.add_argument("--x", type=Path, required=True, help="First variable/data file")
    information_mi.add_argument("--y", type=Path, required=True, help="Second variable/data file")
    information_mi.add_argument("--output", type=Path, default=Path("output/information"), help="Output directory")
    information_profile = information_sub.add_parser("profile", help="Calculate information profile")
    information_profile.add_argument("--sequences", type=Path, required=True, help="Input sequences file")
    information_profile.add_argument("--output", type=Path, default=Path("output/information"), help="Output directory")
    information_profile.add_argument("--k", type=int, default=1, help="K-mer size")
    information_profile.add_argument("--visualize", action="store_true", help="Generate visualization plots")

    # life-events subcommand
    life_events_parser = subparsers.add_parser("life-events", help="Life course analysis operations")
    life_events_sub = life_events_parser.add_subparsers(dest="life_events_cmd")
    life_events_embed = life_events_sub.add_parser("embed", help="Learn event embeddings")
    life_events_embed.add_argument("--input", type=Path, required=True, help="Event sequence file (JSON)")
    life_events_embed.add_argument("--output", type=Path, default=Path("output/life_events/embeddings"), help="Output directory")
    life_events_embed.add_argument("--embedding-dim", type=int, default=100, help="Embedding dimension")
    life_events_embed.add_argument("--window-size", type=int, default=5, help="Context window size")
    life_events_embed.add_argument("--epochs", type=int, default=10, help="Training epochs")
    life_events_predict = life_events_sub.add_parser("predict", help="Predict life outcomes")
    life_events_predict.add_argument("--events", type=Path, required=True, help="Event sequences file")
    life_events_predict.add_argument("--model", type=Path, required=True, help="Pre-trained model path")
    life_events_predict.add_argument("--output", type=Path, default=Path("output/life_events/predictions"), help="Output directory")
    life_events_interpret = life_events_sub.add_parser("interpret", help="Interpret predictions")
    life_events_interpret.add_argument("--model", type=Path, required=True, help="Model path")
    life_events_interpret.add_argument("--sequences", type=Path, required=True, help="Event sequences file")
    life_events_interpret.add_argument("--output", type=Path, default=Path("output/life_events/interpretation"), help="Output directory")

    args = parser.parse_args()

    if args.command == "setup":
        # Execute scripts/setup_uv.sh non-interactively
        root = Path(__file__).resolve().parents[2]
        script = root / "scripts" / "setup_uv.sh"
        cmd = ["bash", str(script)]
        if args.with_amalgkit:
            cmd.append("--with-amalgkit")
        if args.ncbi_email:
            cmd.extend(["--ncbi-email", args.ncbi_email])
        import subprocess

        rc = subprocess.run(cmd).returncode
        sys.exit(rc)

    if args.command == "dna":
        if args.dna_cmd == "fetch":
            # Lazy import here to avoid importing optional Bio dependencies unless needed
            from .dna.genomes import is_valid_assembly_accession

            if not is_valid_assembly_accession(args.assembly):
                print(f"Invalid assembly accession: {args.assembly}")
                sys.exit(2)
            print(f"Validated assembly accession: {args.assembly}")
            # Future: call into actual fetch workflow
            return
        
        if args.dna_cmd == "align":
            from .dna import sequences, alignment
            from .core import io, paths
            
            paths.ensure_directory(args.output)
            
            # Read sequences
            seqs = sequences.read_fasta(args.input)
            seq_list = list(seqs.values())
            
            if args.method == "pairwise":
                if len(seq_list) < 2:
                    print("Error: Need at least 2 sequences for pairwise alignment", file=sys.stderr)
                    sys.exit(1)
                
                # Perform pairwise alignment
                aln_result = alignment.global_pairwise(seq_list[0], seq_list[1]) if args.pairwise_type == "global" else alignment.local_pairwise(seq_list[0], seq_list[1])
                
                result = {
                    "method": args.method,
                    "type": args.pairwise_type,
                    "score": aln_result.get("score", 0),
                    "alignment": aln_result.get("alignment", []),
                }
                
                output_file = args.output / "pairwise_alignment.json"
                io.dump_json(result, output_file)
                print(f"Pairwise alignment complete. Results saved to {output_file}")
                return
            
            elif args.method == "msa":
                from .dna import msa
                
                if len(seq_list) < 2:
                    print("Error: Need at least 2 sequences for MSA", file=sys.stderr)
                    sys.exit(1)
                
                # Perform MSA
                msa_result = msa.multiple_sequence_alignment(seq_list)
                
                result = {
                    "method": "msa",
                    "n_sequences": len(seq_list),
                    "alignment": msa_result,
                }
                
                output_file = args.output / "msa_alignment.json"
                io.dump_json(result, output_file)
                print(f"MSA complete. Results saved to {output_file}")
                return
        
        if args.dna_cmd == "phylogeny":
            from .dna import sequences, phylogeny, alignment
            from .core import io, paths
            
            paths.ensure_directory(args.output)
            
            # Read sequences
            seqs = sequences.read_fasta(args.input)
            seq_list = list(seqs.values())
            seq_ids = list(seqs.keys())
            
            if len(seq_list) < 3:
                print("Error: Need at least 3 sequences for phylogeny", file=sys.stderr)
                sys.exit(1)
            
            # Build distance matrix
            from .dna.distances import pairwise_distances
            dist_matrix = pairwise_distances(seq_list)
            
            # Build tree
            tree = phylogeny.neighbor_joining(dist_matrix, labels=seq_ids)
            
            if args.format == "newick":
                output_file = args.output / "tree.newick"
                output_file.write_text(tree)
                print(f"Phylogenetic tree saved to {output_file}")
            else:
                result = {
                    "method": args.method,
                    "tree": tree,
                    "n_sequences": len(seq_list),
                }
                output_file = args.output / "tree.json"
                io.dump_json(result, output_file)
                print(f"Phylogenetic tree saved to {output_file}")
            return
        
        if args.dna_cmd == "population":
            from .dna import sequences, population
            from .core import io, paths
            
            paths.ensure_directory(args.output)
            
            # Read sequences
            seqs = sequences.read_fasta(args.input)
            seq_list = list(seqs.values())
            
            if len(seq_list) < 2:
                print("Error: Need at least 2 sequences for population analysis", file=sys.stderr)
                sys.exit(1)
            
            results = {}
            
            if "all" in args.statistics or "pi" in args.statistics:
                pi = population.nucleotide_diversity(seq_list)
                results["nucleotide_diversity"] = pi
            
            if "all" in args.statistics or "tajima_d" in args.statistics:
                tajima_d = population.tajimas_d(seq_list)
                results["tajimas_d"] = tajima_d
            
            if "all" in args.statistics or "fst" in args.statistics:
                # For Fst, need two populations - split sequences in half
                mid = len(seq_list) // 2
                pop1 = seq_list[:mid]
                pop2 = seq_list[mid:]
                fst = population.fst(pop1, pop2)
                results["fst"] = fst
            
            output_file = args.output / "population_statistics.json"
            io.dump_json(results, output_file)
            print(f"Population genetics analysis complete. Results saved to {output_file}")
            return
        
        if args.dna_cmd == "variants":
            from .dna import variants
            from .core import io, paths
            
            paths.ensure_directory(args.output)
            
            if args.format == "vcf":
                # Read VCF file
                vcf_data = variants.read_vcf(args.input)
                
                result = {
                    "format": "vcf",
                    "n_variants": len(vcf_data) if isinstance(vcf_data, list) else 0,
                    "variants": vcf_data[:100] if isinstance(vcf_data, list) else [],  # Limit output
                }
                
                output_file = args.output / "variants.json"
                io.dump_json(result, output_file)
                print(f"Variant analysis complete. Results saved to {output_file}")
                return
            else:
                print("Error: BAM format variant calling not yet implemented", file=sys.stderr)
                sys.exit(1)

    if args.command == "rna":
        if args.rna_cmd == "plan":
            cfg = AmalgkitWorkflowConfig(work_dir=Path(args.work_dir), threads=args.threads, species_list=args.species)
            from .rna.workflow import plan_workflow

            steps = plan_workflow(cfg)
            for name, params in steps:
                print(name, params)
            return

        if args.rna_cmd == "run":
            cfg = AmalgkitWorkflowConfig(work_dir=Path(args.work_dir), threads=args.threads, species_list=args.species)
            codes = execute_workflow(cfg, check=args.check)
            print("Return codes:", codes)
            sys.exit(max(codes) if codes else 0)

        if args.rna_cmd == "plan-species":
            base = Path(args.work_dir)
            cfg = AmalgkitWorkflowConfig(work_dir=base, threads=args.threads)
            species = SpeciesProfile(name="", taxon_id=args.taxon_id, tissues=args.tissue or None)
            layout = AmalgkitRunLayout(base_dir=base)
            params_map = build_step_params(species, layout)
            steps = plan_workflow_with_params(cfg, params_map)
            for name, params in steps:
                print(name, params)
            return

        if args.rna_cmd == "plan-config":
            from .rna.workflow import load_workflow_config, plan_workflow

            cfg = load_workflow_config(args.config)
            steps = plan_workflow(cfg)
            for name, params in steps:
                print(name, params)
            return

        if args.rna_cmd == "run-config":
            from .rna.workflow import load_workflow_config

            cfg = load_workflow_config(args.config)
            codes = execute_workflow(cfg, check=args.check)
            print("Return codes:", codes)
            sys.exit(max(codes) if codes else 0)

    if args.command == "protein":
        if args.protein_cmd == "taxon-ids":
            from .protein.proteomes import read_taxon_ids

            ids = read_taxon_ids(Path(args.file))
            print("\n".join(str(i) for i in ids))
            return

        if args.protein_cmd == "comp":
            from .protein.sequences import calculate_aa_composition, parse_fasta

            sequences = parse_fasta(Path(args.fasta))
            for seq_id, seq in sequences.items():
                comp = calculate_aa_composition(seq)
                comp_str = ",".join(f"{aa}:{freq:.3f}" for aa, freq in comp.items() if freq > 0)
                print(f"{seq_id}\t{comp_str}")
            return

        if args.protein_cmd == "rmsd-ca":
            import numpy as np

            from .protein.structure import compute_rmsd_kabsch
            from .protein.structure_io import read_pdb_ca_coordinates

            coords_a = read_pdb_ca_coordinates(Path(args.pdb_a))
            coords_b = read_pdb_ca_coordinates(Path(args.pdb_b))
            if len(coords_a) != len(coords_b):
                print("Error: PDB files have different numbers of CA atoms")
                sys.exit(1)
            rmsd = compute_rmsd_kabsch(np.array(coords_a), np.array(coords_b))
            print(f"{rmsd:.6f}")
            return
        
        if args.protein_cmd == "align":
            from .protein import sequences, alignment
            from .core import io, paths
            
            paths.ensure_directory(args.output)
            
            # Read sequences
            seqs = sequences.parse_fasta(args.input)
            seq_list = list(seqs.values())
            seq_ids = list(seqs.keys())
            
            if len(seq_list) < 2:
                print("Error: Need at least 2 sequences for alignment", file=sys.stderr)
                sys.exit(1)
            
            if args.method == "needleman_wunsch":
                # Pairwise alignment
                score, aln1, aln2 = alignment.needleman_wunsch(seq_list[0], seq_list[1])
                
                result = {
                    "method": args.method,
                    "score": score,
                    "alignment": [aln1, aln2],
                    "sequences": seq_ids[:2],
                }
                
                output_file = args.output / "alignment.json"
                io.dump_json(result, output_file)
                print(f"Alignment complete. Score: {score}, saved to {output_file}")
                return
            else:
                # Pairwise identity for all pairs
                from .protein import pairwise_identity
                
                results = []
                for i, seq1 in enumerate(seq_list):
                    for j, seq2 in enumerate(seq_list[i+1:], i+1):
                        identity = pairwise_identity(seq1, seq2)
                        results.append({
                            "sequence1": seq_ids[i],
                            "sequence2": seq_ids[j],
                            "identity": identity,
                        })
                
                result_df = pd.DataFrame(results)
                output_file = args.output / "pairwise_identity.tsv"
                result_df.to_csv(output_file, sep="\t", index=False)
                print(f"Pairwise identity analysis complete. Saved to {output_file}")
                return
        
        if args.protein_cmd == "structure":
            from .protein import sequences
            from .core import io, paths
            
            paths.ensure_directory(args.output)
            
            # Read sequences
            seqs = sequences.parse_fasta(args.input)
            
            results = {}
            
            for seq_id, seq in seqs.items():
                seq_results = {}
                
                if "all" in args.analyze or "secondary" in args.analyze:
                    try:
                        from .protein import predict_secondary_structure
                        secondary = predict_secondary_structure(seq)
                        seq_results["secondary_structure"] = secondary
                    except ImportError:
                        pass
                
                if "all" in args.analyze or "stability" in args.analyze:
                    try:
                        from .protein import analyze_protein_stability
                        stability = analyze_protein_stability(seq)
                        seq_results["stability"] = stability
                    except ImportError:
                        pass
                
                if "all" in args.analyze or "domains" in args.analyze:
                    try:
                        from .protein import identify_domains
                        domains = identify_domains(seq)
                        seq_results["domains"] = domains
                    except ImportError:
                        pass
                
                if "all" in args.analyze or "ptm" in args.analyze:
                    try:
                        from .protein import analyze_post_translational_modifications
                        ptm = analyze_post_translational_modifications(seq)
                        seq_results["ptm"] = ptm
                    except ImportError:
                        pass
                
                results[seq_id] = seq_results
            
            output_file = args.output / "structure_analysis.json"
            io.dump_json(results, output_file, indent=2)
            print(f"Structure analysis complete. Saved to {output_file}")
            return
        
        if args.protein_cmd == "domains":
            from .protein import sequences
            from .core import io, paths
            
            paths.ensure_directory(args.output)
            
            # Read sequences
            seqs = sequences.parse_fasta(args.input)
            
            domain_results = []
            
            for seq_id, seq in seqs.items():
                if args.uniprot:
                    # Try to fetch from UniProt (would need sequence ID mapping)
                    try:
                        from .protein import fetch_interpro_domains
                        # Would need UniProt accession - placeholder
                        domains = []
                    except ImportError:
                        domains = []
                else:
                    try:
                        from .protein import identify_domains
                        domains = identify_domains(seq)
                    except ImportError:
                        domains = []
                
                domain_results.append({
                    "sequence_id": seq_id,
                    "n_domains": len(domains),
                    "domains": domains,
                })
            
            result_df = pd.DataFrame(domain_results)
            output_file = args.output / "domains.json"
            io.dump_json(domain_results, output_file, indent=2)
            print(f"Domain annotation complete. Saved to {output_file}")
            return
        
        if args.protein_cmd == "families":
            from .protein import sequences
            from .core import io, paths
            
            paths.ensure_directory(args.output)
            
            # Read sequences
            seqs = sequences.parse_fasta(args.input)
            
            family_results = []
            
            for seq_id, seq in seqs.items():
                try:
                    from .protein import predict_protein_family
                    families = predict_protein_family(seq)
                    family_results.append({
                        "sequence_id": seq_id,
                        "predicted_families": families[:5],  # Top 5
                    })
                except ImportError:
                    family_results.append({
                        "sequence_id": seq_id,
                        "predicted_families": [],
                    })
            
            output_file = args.output / "families.json"
            io.dump_json(family_results, output_file, indent=2)
            print(f"Family classification complete. Saved to {output_file}")
            return

    if args.command == "math":
        if args.math_cmd == "popgen":
            from .dna import sequences
            from .core import io, paths
            
            paths.ensure_directory(args.output)
            
            if args.input:
                seqs = sequences.read_fasta(args.input)
                seq_list = list(seqs.values())
            else:
                # Generate synthetic data
                from .simulation.popgen import generate_population_sequences
                seq_list = generate_population_sequences(n_sequences=10, sequence_length=1000)
            
            results = {}
            
            if "all" in args.analysis or "tajima_d" in args.analysis:
                from .dna.population import tajimas_d
                tajima_d = tajimas_d(seq_list)
                results["tajimas_d"] = tajima_d
            
            if "all" in args.analysis or "fst" in args.analysis:
                from .dna.population import fst
                if len(seq_list) >= 4:
                    mid = len(seq_list) // 2
                    pop1 = seq_list[:mid]
                    pop2 = seq_list[mid:]
                    fst_value = fst(pop1, pop2)
                    results["fst"] = fst_value
            
            if "all" in args.analysis or "ld" in args.analysis:
                # LD calculation requires genotype data (not just sequences)
                # Would need phased haplotypes or genotype calls
                results["ld_analysis"] = {
                    "note": "LD calculation requires genotype data (phased haplotypes or VCF)",
                    "available_functions": "Use metainformant.math.ld.ld_coefficients() with haplotype frequencies",
                }
            
            if "all" in args.analysis or "coalescent" in args.analysis:
                from .math.coalescent import expected_time_to_mrca, expected_total_branch_length
                # Calculate expected coalescent statistics
                Ne = 1000.0  # Default effective population size
                t_mrca = expected_time_to_mrca(len(seq_list), Ne)
                total_length = expected_total_branch_length(len(seq_list), Ne)
                results["coalescent"] = {
                    "expected_time_to_mrca": float(t_mrca),
                    "expected_total_branch_length": float(total_length),
                    "effective_population_size": Ne,
                }
            
            output_file = args.output / "popgen_analysis.json"
            io.dump_json(results, output_file)
            print(f"Population genetics analysis complete. Saved to {output_file}")
            return
        
        if args.math_cmd == "coalescent":
            from .math.coalescent import (
                expected_coalescent_waiting_times,
                expected_time_to_mrca,
                expected_total_branch_length,
            )
            from .math.demography import (
                bottleneck_effective_size,
                exponential_growth_effective_size,
            )
            from .core import io, paths
            
            paths.ensure_directory(args.output)
            
            # Calculate effective population size based on model
            if args.model == "constant":
                Ne = 1000.0
            elif args.model == "exponential":
                Ne = exponential_growth_effective_size(
                    current_size=1000,
                    growth_rate=0.01,
                )
            else:  # bottleneck
                Ne = bottleneck_effective_size(
                    pre_bottleneck_size=1000,
                    bottleneck_size=100,
                    bottleneck_duration=0.1,
                )
            
            # Calculate coalescent statistics
            t_mrca = expected_time_to_mrca(args.n_samples, Ne)
            total_length = expected_total_branch_length(args.n_samples, Ne)
            waiting_times = expected_coalescent_waiting_times(args.n_samples, Ne)
            
            results = {
                "model": args.model,
                "n_samples": args.n_samples,
                "n_loci": args.n_loci,
                "effective_population_size": float(Ne),
                "expected_time_to_mrca": float(t_mrca),
                "expected_total_branch_length": float(total_length),
                "expected_waiting_times": [float(t) for t in waiting_times],
            }
            
            output_file = args.output / "coalescent_results.json"
            io.dump_json(results, output_file)
            print(f"Coalescent analysis complete. Saved to {output_file}")
            return
        
        # Handle selection experiments (existing)
        func = getattr(args, "func", None)
        if callable(func):
            func(args)
            return

    if args.command == "gwas":
        if args.gwas_cmd == "run":
            from .gwas.workflow import execute_gwas_workflow
            from .gwas.config import load_gwas_config

            cfg = load_gwas_config(args.config)
            results = execute_gwas_workflow(cfg, check=args.check)
            if results.get("status") == "completed":
                print(f"GWAS workflow completed successfully. Results in {cfg.work_dir}")
                sys.exit(0)
            elif results.get("status") == "failed":
                print(f"GWAS workflow failed: {results.get('error', 'Unknown error')}")
                sys.exit(1)
            else:
                print(f"GWAS workflow status: {results.get('status', 'unknown')}")
                sys.exit(0)

    if args.command == "ontology":
        if args.ontology_cmd == "run":
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "ontology" / "run_ontology_analysis.py"
            cmd = ["python3", str(script), f"--output={args.output}"]
            if args.go:
                cmd.append(f"--go={args.go}")
            if args.query_term:
                cmd.extend(["--query-term", args.query_term])
            if args.ancestors:
                cmd.append("--ancestors")
            if args.descendants:
                cmd.append("--descendants")
            sys.exit(subprocess.run(cmd).returncode)

    if args.command == "phenotype":
        if args.phenotype_cmd == "run":
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "phenotype" / "run_phenotype_analysis.py"
            cmd = ["python3", str(script), f"--input={args.input}", f"--output={args.output}"]
            if args.analyze_statistics:
                cmd.append("--analyze-statistics")
            if args.analyze_correlations:
                cmd.append("--analyze-correlations")
            sys.exit(subprocess.run(cmd).returncode)

    if args.command == "networks":
        if args.networks_cmd == "run":
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "networks" / "run_network_analysis.py"
            cmd = ["python3", str(script), f"--input={args.input}", f"--output={args.output}"]
            if args.analyze_metrics:
                cmd.append("--analyze-metrics")
            if args.detect_communities:
                cmd.append("--detect-communities")
            if args.analyze_centrality:
                cmd.append("--analyze-centrality")
            sys.exit(subprocess.run(cmd).returncode)

    if args.command == "multiomics":
        if args.multiomics_cmd == "run":
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "multiomics" / "run_multiomics_integration.py"
            cmd = ["python3", str(script), f"--output={args.output}"]
            if args.genomics:
                cmd.append(f"--genomics={args.genomics}")
            if args.transcriptomics:
                cmd.append(f"--transcriptomics={args.transcriptomics}")
            if args.proteomics:
                cmd.append(f"--proteomics={args.proteomics}")
            if args.joint_pca:
                cmd.append("--joint-pca")
            if args.joint_nmf:
                cmd.append("--joint-nmf")
            if args.canonical_correlation:
                cmd.append("--canonical-correlation")
            sys.exit(subprocess.run(cmd).returncode)

    if args.command == "singlecell":
        if args.singlecell_cmd == "run":
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "singlecell" / "run_singlecell_analysis.py"
            cmd = ["python3", str(script), f"--input={args.input}", f"--output={args.output}"]
            if args.qc:
                cmd.append("--qc")
            if args.normalize:
                cmd.append("--normalize")
            if args.cluster:
                cmd.append("--cluster")
            sys.exit(subprocess.run(cmd).returncode)

    if args.command == "quality":
        if args.quality_cmd == "run":
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "quality" / "run_quality_control.py"
            cmd = ["python3", str(script), f"--output={args.output}"]
            if args.fastq:
                cmd.append(f"--fastq={args.fastq}")
            if args.analyze_fastq:
                cmd.append("--analyze-fastq")
            if args.detect_contamination:
                cmd.append("--detect-contamination")
            sys.exit(subprocess.run(cmd).returncode)

    if args.command == "simulation":
        if args.simulation_cmd == "run":
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "simulation" / "run_simulation.py"
            cmd = ["python3", str(script), f"--model={args.model}", f"--output={args.output}"]
            if args.n:
                cmd.append(f"--n={args.n}")
            if args.steps:
                cmd.append(f"--steps={args.steps}")
            sys.exit(subprocess.run(cmd).returncode)

    if args.command == "visualization":
        if args.visualization_cmd == "run":
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "visualization" / "run_visualization.py"
            cmd = ["python3", str(script), f"--input={args.input}", f"--plot-type={args.plot_type}", f"--output={args.output}"]
            sys.exit(subprocess.run(cmd).returncode)

    if args.command == "epigenome":
        if args.epigenome_cmd == "run":
            from .epigenome import run_methylation_workflow, run_chipseq_workflow, run_atacseq_workflow
            from .core import paths
            
            paths.ensure_directory(args.output)
            
            if args.methylation:
                # Run methylation workflow
                try:
                    results = run_methylation_workflow(
                        args.methylation,
                        args.output,
                        compute_beta=args.compute_beta,
                    )
                    print(f"Methylation analysis complete. Results in {args.output}")
                    return
                except Exception as e:
                    print(f"Error running methylation workflow: {e}", file=sys.stderr)
                    sys.exit(1)
            
            if args.bedgraph:
                # Determine if ChIP-seq or ATAC-seq based on file or user input
                # For now, try ChIP-seq
                from .epigenome.tracks import read_bedgraph
                try:
                    signal = read_bedgraph(args.bedgraph)
                    
                    if not signal.empty:
                        # Get chromosome from data
                        chroms = signal["chrom"].unique()
                        if len(chroms) > 0:
                            chrom = chroms[0]
                            results = run_chipseq_workflow(args.bedgraph, args.output, chrom=chrom)
                            print(f"ChIP-seq analysis complete. Results in {args.output}")
                            return
                    else:
                        print("Warning: bedGraph file is empty", file=sys.stderr)
                except Exception as e:
                    print(f"Error reading bedGraph file: {e}", file=sys.stderr)
                    sys.exit(1)
            
            print("Error: Specify --methylation or --bedgraph", file=sys.stderr)
            sys.exit(1)

    if args.command == "ecology":
        if args.ecology_cmd == "run":
            from .ecology import run_community_analysis_workflow
            from .core import paths
            
            paths.ensure_directory(args.output)
            
            try:
                results = run_community_analysis_workflow(
                    args.input,
                    args.output,
                    calculate_diversity=args.diversity,
                    calculate_beta=args.beta_diversity,
                )
                print(f"Ecology analysis complete. Results in {args.output}")
                return
            except Exception as e:
                print(f"Error running ecology workflow: {e}", file=sys.stderr)
                sys.exit(1)

    if args.command == "ml":
        if args.ml_cmd == "run":
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "ml" / "run_ml_pipeline.py"
            cmd = ["python3", str(script), f"--features={args.features}", f"--output={args.output}"]
            if args.labels:
                cmd.append(f"--labels={args.labels}")
            if args.classify:
                cmd.append("--classify")
            if args.regress:
                cmd.append("--regress")
            if args.feature_selection:
                cmd.append("--feature-selection")
            sys.exit(subprocess.run(cmd).returncode)

    if args.command == "information":
        from metainformant.core import io, paths
        from metainformant.information import (
            analyze_sequence_information,
            information_profile,
            mutual_information,
            shannon_entropy,
        )
        from metainformant.information.visualization import plot_information_profile

        if args.information_cmd == "entropy":
            # Read input data
            input_path = Path(args.input)
            output_dir = Path(args.output)
            paths.ensure_directory(output_dir)

            if input_path.suffix == ".fasta" or input_path.suffix == ".fa":
                # Read FASTA file
                from metainformant.dna import sequences

                seqs = sequences.read_fasta(input_path)
                sequences_list = list(seqs.values())
            else:
                # Assume text file with one sequence per line
                with open(input_path) as f:
                    sequences_list = [line.strip() for line in f if line.strip()]

            # Calculate entropy for each sequence
            results = []
            for seq in sequences_list:
                analysis = analyze_sequence_information(seq, k_values=[args.k])
                entropy = analysis["kmer_analyses"][args.k]["entropy"]
                results.append({"sequence": seq[:50] + "..." if len(seq) > 50 else seq, "entropy": entropy})

            # Save results
            output_file = output_dir / "entropy_results.json"
            io.dump_json(results, output_file)
            print(f"Entropy analysis complete. Results saved to {output_file}")
            return

        if args.information_cmd == "mutual-information":
            # Read input data
            x_path = Path(args.x)
            y_path = Path(args.y)
            output_dir = Path(args.output)
            paths.ensure_directory(output_dir)

            # Read sequences or data
            with open(x_path) as f:
                x_data = [line.strip() for line in f if line.strip()]
            with open(y_path) as f:
                y_data = [line.strip() for line in f if line.strip()]

            if len(x_data) != len(y_data):
                print(f"Error: X and Y must have the same length ({len(x_data)} vs {len(y_data)})")
                sys.exit(1)

            # Calculate mutual information
            mi = mutual_information(x_data, y_data)

            # Save results
            results = {"mutual_information": mi, "n_samples": len(x_data)}
            output_file = output_dir / "mutual_information_results.json"
            io.dump_json(results, output_file)
            print(f"Mutual information: {mi:.6f} bits")
            print(f"Results saved to {output_file}")
            return

        if args.information_cmd == "profile":
            # Read sequences
            sequences_path = Path(args.sequences)
            output_dir = Path(args.output)
            paths.ensure_directory(output_dir)

            if sequences_path.suffix == ".fasta" or sequences_path.suffix == ".fa":
                from metainformant.dna import sequences

                seqs = sequences.read_fasta(sequences_path)
                sequences_list = list(seqs.values())
            else:
                with open(sequences_path) as f:
                    sequences_list = [line.strip() for line in f if line.strip()]

            # Calculate information profile
            profile = information_profile(sequences_list, k=args.k)

            # Save results
            output_file = output_dir / "information_profile.json"
            io.dump_json(profile, output_file)
            print(f"Information profile complete. Entropy: {profile['entropy']:.6f} bits")
            print(f"Results saved to {output_file}")

            # Generate visualization if requested
            if args.visualize:
                try:
                    plot_info = plot_information_profile(profile, output_path=output_dir / "profile.png")
                    print(f"Visualization saved to {plot_info['output_path']}")
                except Exception as e:
                    print(f"Warning: Could not generate visualization: {e}")

            return

    if args.command == "life-events":
        if args.life_events_cmd == "embed":
            from .life_events import EventDatabase, learn_event_embeddings, load_sequences_from_json
            from .core import io

            # Validate input file exists
            if not args.input.exists():
                print(f"Error: Input file not found: {args.input}", file=sys.stderr)
                sys.exit(1)

            try:
                # Load event sequences using utility function
                sequences = load_sequences_from_json(args.input)
                database = EventDatabase(sequences=sequences)
            except Exception as e:
                print(f"Error loading event sequences: {e}", file=sys.stderr)
                sys.exit(1)

            if not database.sequences:
                print("Error: No sequences found in input file", file=sys.stderr)
                sys.exit(1)

            # Convert to token sequences
            sequences_tokens = []
            for seq in database.sequences:
                tokens = [f"{e.domain}:{e.event_type}" for e in seq.events]
                sequences_tokens.append(tokens)

            if not sequences_tokens or all(not tokens for tokens in sequences_tokens):
                print("Error: No valid event sequences found", file=sys.stderr)
                sys.exit(1)

            try:
                # Learn embeddings
                embeddings = learn_event_embeddings(
                    sequences_tokens,
                    embedding_dim=args.embedding_dim,
                    window_size=args.window_size,
                    epochs=args.epochs
                )

                # Save embeddings
                from .core.paths import prepare_file_path
                output_file = args.output / "embeddings.json"
                prepare_file_path(output_file)
                embeddings_dict = {k: v.tolist() for k, v in embeddings.items()}
                io.dump_json(embeddings_dict, output_file, indent=2)
                print(f"Learned {len(embeddings)} event embeddings. Saved to {output_file}")
                return
            except Exception as e:
                print(f"Error learning embeddings: {e}", file=sys.stderr)
                sys.exit(1)

        if args.life_events_cmd == "predict":
            import numpy as np
            from .life_events import EventDatabase, EventSequencePredictor, load_sequences_from_json
            from .core import io, paths

            # Validate input file exists
            if not args.events.exists():
                print(f"Error: Events file not found: {args.events}", file=sys.stderr)
                sys.exit(1)

            # Validate model file exists
            if args.model is None or not args.model.exists():
                print(f"Error: Model file not found: {args.model}", file=sys.stderr)
                print("Note: Use 'life-events embed' to train a model first, or provide --model path", file=sys.stderr)
                sys.exit(1)

            try:
                # Load event sequences
                sequences = load_sequences_from_json(args.events)
                database = EventDatabase(sequences=sequences)
            except Exception as e:
                print(f"Error loading event sequences: {e}", file=sys.stderr)
                sys.exit(1)

            if not database.sequences:
                print("Error: No sequences found in events file", file=sys.stderr)
                sys.exit(1)

            try:
                # Load trained model
                predictor = EventSequencePredictor.load_model(args.model)
            except Exception as e:
                print(f"Error loading model: {e}", file=sys.stderr)
                sys.exit(1)

            # Convert to token sequences
            sequences_tokens = []
            for seq in database.sequences:
                tokens = [f"{e.domain}:{e.event_type}" for e in seq.events]
                sequences_tokens.append(tokens)

            try:
                # Make predictions
                predictions = predictor.predict(sequences_tokens)
                
                # Prepare output
                output_dir = Path(args.output)
                paths.ensure_directory(output_dir)
                
                # Create predictions dictionary with sequence IDs
                predictions_dict = {
                    "n_sequences": len(sequences),
                    "model_path": str(args.model),
                    "model_type": predictor.model_type,
                    "task_type": predictor.task_type,
                    "predictions": []
                }
                
                # Add predictions with sequence IDs
                for i, seq in enumerate(database.sequences):
                    pred_entry = {
                        "sequence_id": seq.person_id,
                        "prediction": float(predictions[i])
                    }
                    
                    # Add probabilities if classification task
                    if predictor.task_type == "classification":
                        try:
                            proba = predictor.predict_proba([sequences_tokens[i]])
                            pred_entry["probabilities"] = {
                                str(cls): float(prob) 
                                for cls, prob in zip(predictor.classes_, proba[0])
                            }
                        except Exception:
                            pass  # Skip if predict_proba fails
                    
                    predictions_dict["predictions"].append(pred_entry)
                
                # Save predictions
                predictions_file = output_dir / "predictions.json"
                io.dump_json(predictions_dict, predictions_file, indent=2)
                
                print(f"Made predictions for {len(sequences)} sequences")
                print(f"Predictions saved to {predictions_file}")
                
                # Print summary
                if predictor.task_type == "classification":
                    unique_preds, counts = np.unique(predictions, return_counts=True)
                    print("\nPrediction summary:")
                    for pred, count in zip(unique_preds, counts):
                        print(f"  Class {pred}: {count} sequences")
                else:
                    print(f"\nPrediction statistics:")
                    print(f"  Mean: {np.mean(predictions):.4f}")
                    print(f"  Std: {np.std(predictions):.4f}")
                    print(f"  Min: {np.min(predictions):.4f}")
                    print(f"  Max: {np.max(predictions):.4f}")
                
                return
            except Exception as e:
                print(f"Error making predictions: {e}", file=sys.stderr)
                import traceback
                traceback.print_exc()
                sys.exit(1)

        if args.life_events_cmd == "interpret":
            from .life_events import (
                EventDatabase,
                EventSequencePredictor,
                event_importance,
                feature_attribution,
                load_sequences_from_json,
                temporal_patterns,
            )
            from .core import io, paths

            # Validate model file exists
            if not args.model.exists():
                print(f"Error: Model file not found: {args.model}", file=sys.stderr)
                sys.exit(1)

            # Validate sequences file exists
            if not args.sequences.exists():
                print(f"Error: Sequences file not found: {args.sequences}", file=sys.stderr)
                sys.exit(1)

            try:
                # Load model
                predictor = EventSequencePredictor.load_model(args.model)
            except Exception as e:
                print(f"Error loading model: {e}", file=sys.stderr)
                sys.exit(1)

            try:
                # Load sequences
                sequences = load_sequences_from_json(args.sequences)
                database = EventDatabase(sequences=sequences)
            except Exception as e:
                print(f"Error loading sequences: {e}", file=sys.stderr)
                sys.exit(1)

            if not database.sequences:
                print("Error: No sequences found in sequences file", file=sys.stderr)
                sys.exit(1)

            # Convert to token sequences
            sequences_tokens = []
            for seq in database.sequences:
                tokens = [f"{e.domain}:{e.event_type}" for e in seq.events]
                sequences_tokens.append(tokens)

            # Get event embeddings from model
            if not hasattr(predictor, "event_embeddings") or not predictor.event_embeddings:
                print("Error: Model does not contain event embeddings required for interpretation", file=sys.stderr)
                sys.exit(1)

            event_embeddings = predictor.event_embeddings

            try:
                # Prepare output directory
                output_dir = Path(args.output)
                paths.ensure_directory(output_dir)

                # Compute interpretations
                print("Computing event importance...")
                importance = event_importance(predictor, sequences_tokens, event_embeddings, method="permutation")

                # Get predictions for temporal patterns analysis
                print("Making predictions for temporal analysis...")
                predictions = predictor.predict(sequences_tokens)

                print("Computing temporal patterns...")
                temporal = temporal_patterns(sequences_tokens, predictions)

                print("Computing feature attribution...")
                attribution = feature_attribution(predictor, sequences_tokens, event_embeddings, use_shap=False)

                # Create interpretation report
                interpretation_report = {
                    "model_path": str(args.model),
                    "model_type": predictor.model_type,
                    "task_type": predictor.task_type,
                    "n_sequences": len(sequences_tokens),
                    "interpretations": {
                        "event_importance": importance,
                        "temporal_patterns": temporal,
                        "feature_attribution": attribution,
                    }
                }

                # Save report
                report_file = output_dir / "interpretation_report.json"
                io.dump_json(interpretation_report, report_file, indent=2)

                print(f"\nInterpretation complete. Report saved to {report_file}")

                # Print top important events
                sorted_importance = sorted(importance.items(), key=lambda x: x[1], reverse=True)
                print("\nTop 10 most important events:")
                for i, (event, score) in enumerate(sorted_importance[:10], 1):
                    print(f"  {i}. {event}: {score:.4f}")

                # Try to create visualization if matplotlib available
                try:
                    from .life_events import plot_prediction_importance
                    if plot_prediction_importance is not None:
                        viz_file = output_dir / "importance_plot.png"
                        plot_prediction_importance(importance, top_n=20, output_path=viz_file)
                        print(f"\nVisualization saved to {viz_file}")
                except Exception as e:
                    # Visualization is optional
                    pass

                return
            except Exception as e:
                print(f"Error computing interpretations: {e}", file=sys.stderr)
                import traceback
                traceback.print_exc()
                sys.exit(1)

    if args.command == "tests":
        exit_code = run_all_tests(args.pytest_args)
        sys.exit(exit_code)

    parser.print_help()


if __name__ == "__main__":
    main()
