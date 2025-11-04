from __future__ import annotations

import argparse
import sys
from pathlib import Path

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

    # math subcommand
    math_parser = subparsers.add_parser("math", help="Math experiments and demos")
    math_sub = math_parser.add_subparsers(dest="math_cmd")
    # selection sub-tree (wire lazily to keep import cost low)
    from .math.selection_experiments.cli import add_math_selection_subparser as _add_sel

    _add_sel(math_sub)

    # gwas subcommand
    gwas_parser = subparsers.add_parser("gwas", help="GWAS (Genome-Wide Association Studies) workflow")
    gwas_sub = gwas_parser.add_subparsers(dest="gwas_cmd")
    gwas_run = gwas_sub.add_parser("run", help="Run GWAS workflow from configuration file")
    gwas_run.add_argument("--config", required=True, help="Path to GWAS YAML/TOML/JSON config file")
    gwas_run.add_argument("--check", action="store_true", help="Validate configuration only, do not execute")

    # ontology subcommand
    ontology_parser = subparsers.add_parser("ontology", help="Ontology analysis operations")
    ontology_sub = ontology_parser.add_subparsers(dest="ontology_cmd")
    ontology_run = ontology_sub.add_parser("run", help="Run ontology analysis workflow")
    ontology_run.add_argument("--go", type=Path, help="Gene Ontology OBO file")
    ontology_run.add_argument("--output", type=Path, default=Path("output/ontology"), help="Output directory")
    ontology_run.add_argument("--query-term", type=str, help="Query specific GO term")
    ontology_run.add_argument("--ancestors", action="store_true", help="Get ancestors")
    ontology_run.add_argument("--descendants", action="store_true", help="Get descendants")

    # phenotype subcommand
    phenotype_parser = subparsers.add_parser("phenotype", help="Phenotype analysis operations")
    phenotype_sub = phenotype_parser.add_subparsers(dest="phenotype_cmd")
    phenotype_run = phenotype_sub.add_parser("run", help="Run phenotype analysis workflow")
    phenotype_run.add_argument("--input", type=Path, required=True, help="Input phenotype data file")
    phenotype_run.add_argument("--output", type=Path, default=Path("output/phenotype"), help="Output directory")
    phenotype_run.add_argument("--analyze-statistics", action="store_true", help="Calculate statistics")
    phenotype_run.add_argument("--analyze-correlations", action="store_true", help="Calculate correlations")

    # networks subcommand
    networks_parser = subparsers.add_parser("networks", help="Network analysis operations")
    networks_sub = networks_parser.add_subparsers(dest="networks_cmd")
    networks_run = networks_sub.add_parser("run", help="Run network analysis workflow")
    networks_run.add_argument("--input", type=Path, required=True, help="Input interaction file")
    networks_run.add_argument("--output", type=Path, default=Path("output/networks"), help="Output directory")
    networks_run.add_argument("--analyze-metrics", action="store_true", help="Calculate network metrics")
    networks_run.add_argument("--detect-communities", action="store_true", help="Detect communities")
    networks_run.add_argument("--analyze-centrality", action="store_true", help="Calculate centrality")

    # multiomics subcommand
    multiomics_parser = subparsers.add_parser("multiomics", help="Multi-omics integration operations")
    multiomics_sub = multiomics_parser.add_subparsers(dest="multiomics_cmd")
    multiomics_run = multiomics_sub.add_parser("run", help="Run multi-omics integration workflow")
    multiomics_run.add_argument("--genomics", type=Path, help="Genomics data file")
    multiomics_run.add_argument("--transcriptomics", type=Path, help="Transcriptomics data file")
    multiomics_run.add_argument("--proteomics", type=Path, help="Proteomics data file")
    multiomics_run.add_argument("--output", type=Path, default=Path("output/multiomics"), help="Output directory")
    multiomics_run.add_argument("--joint-pca", action="store_true", help="Perform joint PCA")
    multiomics_run.add_argument("--joint-nmf", action="store_true", help="Perform joint NMF")
    multiomics_run.add_argument("--canonical-correlation", action="store_true", help="Perform CCA")

    # singlecell subcommand
    singlecell_parser = subparsers.add_parser("singlecell", help="Single-cell analysis operations")
    singlecell_sub = singlecell_parser.add_subparsers(dest="singlecell_cmd")
    singlecell_run = singlecell_sub.add_parser("run", help="Run single-cell analysis workflow")
    singlecell_run.add_argument("--input", type=Path, required=True, help="Input count matrix file")
    singlecell_run.add_argument("--output", type=Path, default=Path("output/singlecell"), help="Output directory")
    singlecell_run.add_argument("--qc", action="store_true", help="Perform quality control")
    singlecell_run.add_argument("--normalize", action="store_true", help="Normalize counts")
    singlecell_run.add_argument("--cluster", action="store_true", help="Perform clustering")

    # quality subcommand
    quality_parser = subparsers.add_parser("quality", help="Quality control operations")
    quality_sub = quality_parser.add_subparsers(dest="quality_cmd")
    quality_run = quality_sub.add_parser("run", help="Run quality control workflow")
    quality_run.add_argument("--fastq", type=Path, help="Input FASTQ file")
    quality_run.add_argument("--output", type=Path, default=Path("output/quality"), help="Output directory")
    quality_run.add_argument("--analyze-fastq", action="store_true", help="Analyze FASTQ quality")
    quality_run.add_argument("--detect-contamination", action="store_true", help="Detect contamination")

    # simulation subcommand
    simulation_parser = subparsers.add_parser("simulation", help="Simulation operations")
    simulation_sub = simulation_parser.add_subparsers(dest="simulation_cmd")
    simulation_run = simulation_sub.add_parser("run", help="Run simulation workflow")
    simulation_run.add_argument("--model", type=str, required=True, choices=["sequences", "agents", "expression"], help="Simulation model")
    simulation_run.add_argument("--output", type=Path, default=Path("output/simulation"), help="Output directory")
    simulation_run.add_argument("--n", type=int, default=100, help="Number of sequences (sequences model)")
    simulation_run.add_argument("--steps", type=int, default=100, help="Number of steps (agents model)")

    # visualization subcommand
    visualization_parser = subparsers.add_parser("visualization", help="Visualization operations")
    visualization_sub = visualization_parser.add_subparsers(dest="visualization_cmd")
    visualization_run = visualization_sub.add_parser("run", help="Run visualization workflow")
    visualization_run.add_argument("--input", type=Path, required=True, help="Input data file")
    visualization_run.add_argument("--plot-type", type=str, required=True, choices=["lineplot", "heatmap", "animation", "histogram"], help="Plot type")
    visualization_run.add_argument("--output", type=Path, default=Path("output/visualization"), help="Output directory")

    # epigenome subcommand
    epigenome_parser = subparsers.add_parser("epigenome", help="Epigenome analysis operations")
    epigenome_sub = epigenome_parser.add_subparsers(dest="epigenome_cmd")
    epigenome_run = epigenome_sub.add_parser("run", help="Run epigenome analysis workflow")
    epigenome_run.add_argument("--methylation", type=Path, help="Methylation data file")
    epigenome_run.add_argument("--bedgraph", type=Path, help="BedGraph track file")
    epigenome_run.add_argument("--output", type=Path, default=Path("output/epigenome"), help="Output directory")
    epigenome_run.add_argument("--compute-beta", action="store_true", help="Compute beta values")

    # ecology subcommand
    ecology_parser = subparsers.add_parser("ecology", help="Ecology analysis operations")
    ecology_sub = ecology_parser.add_subparsers(dest="ecology_cmd")
    ecology_run = ecology_sub.add_parser("run", help="Run ecology analysis workflow")
    ecology_run.add_argument("--input", type=Path, required=True, help="Input abundance table")
    ecology_run.add_argument("--output", type=Path, default=Path("output/ecology"), help="Output directory")
    ecology_run.add_argument("--diversity", action="store_true", help="Calculate diversity indices")
    ecology_run.add_argument("--beta-diversity", action="store_true", help="Calculate beta diversity")

    # ml subcommand
    ml_parser = subparsers.add_parser("ml", help="Machine learning operations")
    ml_sub = ml_parser.add_subparsers(dest="ml_cmd")
    ml_run = ml_sub.add_parser("run", help="Run ML pipeline workflow")
    ml_run.add_argument("--features", type=Path, required=True, help="Feature matrix file")
    ml_run.add_argument("--labels", type=Path, help="Labels file")
    ml_run.add_argument("--output", type=Path, default=Path("output/ml"), help="Output directory")
    ml_run.add_argument("--classify", action="store_true", help="Perform classification")
    ml_run.add_argument("--regress", action="store_true", help="Perform regression")
    ml_run.add_argument("--feature-selection", action="store_true", help="Perform feature selection")

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
    life_events_predict.add_argument("--model", type=Path, help="Pre-trained model path")
    life_events_predict.add_argument("--outcome", choices=["mortality", "personality", "health"], help="Outcome type to predict")
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

    if args.command == "dna" and args.dna_cmd == "fetch":
        # Lazy import here to avoid importing optional Bio dependencies unless needed
        from .dna.genomes import is_valid_assembly_accession

        if not is_valid_assembly_accession(args.assembly):
            print(f"Invalid assembly accession: {args.assembly}")
            sys.exit(2)
        print(f"Validated assembly accession: {args.assembly}")
        # Future: call into actual fetch workflow
        return

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

    if args.command == "math":
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
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "epigenome" / "run_epigenome_analysis.py"
            cmd = ["python3", str(script), f"--output={args.output}"]
            if args.methylation:
                cmd.append(f"--methylation={args.methylation}")
            if args.bedgraph:
                cmd.append(f"--bedgraph={args.bedgraph}")
            if args.compute_beta:
                cmd.append("--compute-beta")
            sys.exit(subprocess.run(cmd).returncode)

    if args.command == "ecology":
        if args.ecology_cmd == "run":
            import subprocess
            root = Path(__file__).resolve().parents[2]
            script = root / "scripts" / "ecology" / "run_ecology_analysis.py"
            cmd = ["python3", str(script), f"--input={args.input}", f"--output={args.output}"]
            if args.diversity:
                cmd.append("--diversity")
            if args.beta_diversity:
                cmd.append("--beta-diversity")
            sys.exit(subprocess.run(cmd).returncode)

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
            from .life_events import EventDatabase, learn_event_embeddings
            from .core import io, paths

            # Load event sequences
            data = io.load_json(args.input)
            if isinstance(data, dict) and "sequences" in data:
                database = EventDatabase.from_dict(data)
            else:
                # Assume it's a list of sequences
                from .life_events.events import EventSequence
                sequences = [EventSequence.from_dict(s) if isinstance(s, dict) else s for s in data]
                database = EventDatabase(sequences=sequences)

            # Convert to token sequences
            sequences_tokens = []
            for seq in database.sequences:
                tokens = [f"{e.domain}:{e.event_type}" for e in seq.events]
                sequences_tokens.append(tokens)

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
            print(f"Learned {len(embeddings)} event embeddings. Saved to {args.output}")
            return

        if args.life_events_cmd == "predict":
            import numpy as np
            from .life_events import EventDatabase, EventSequencePredictor
            from .core import io

            # Load event sequences
            data = io.load_json(args.events)
            if isinstance(data, dict) and "sequences" in data:
                database = EventDatabase.from_dict(data)
            else:
                from .life_events.events import EventSequence
                sequences = [EventSequence.from_dict(s) if isinstance(s, dict) else s for s in data]
                database = EventDatabase(sequences=sequences)

            # Convert to token sequences
            sequences_tokens = []
            for seq in database.sequences:
                tokens = [f"{e.domain}:{e.event_type}" for e in seq.events]
                sequences_tokens.append(tokens)

            # Note: Full model persistence would require saving/loading model weights
            # For now, demonstrate interface
            print("Note: Prediction requires a trained model.")
            print("Use the analyze_life_course workflow function to train models.")
            print("Full model persistence will be implemented in future updates.")
            sys.exit(0)

        if args.life_events_cmd == "interpret":
            print("Interpretation functionality will be implemented in Phase 3")
            sys.exit(0)

    if args.command == "tests":
        exit_code = run_all_tests(args.pytest_args)
        sys.exit(exit_code)

    parser.print_help()


if __name__ == "__main__":
    main()
