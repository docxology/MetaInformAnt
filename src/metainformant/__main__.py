"""CLI entry point for METAINFORMANT.

This module provides the command-line interface for the METAINFORMANT
bioinformatics toolkit. Run with --help for usage information.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

from . import __version__


def main() -> int:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="metainformant",
        description="METAINFORMANT: Comprehensive Bioinformatics Toolkit for Multi-Omic Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  metainformant --version
  metainformant --help

For detailed usage of specific modules, import them directly in Python:
  python -c "from metainformant.dna.sequence import composition; help(composition.gc_content)"
        """,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"METAINFORMANT {__version__}",
        help="Show version information and exit",
    )

    parser.add_argument(
        "--modules",
        action="store_true",
        help="List available modules",
    )

    subparsers = parser.add_subparsers(dest="command")

    # Protein subcommands
    protein_parser = subparsers.add_parser("protein", help="Protein analysis commands")
    protein_sub = protein_parser.add_subparsers(dest="protein_command")

    # protein taxon-ids
    taxon_parser = protein_sub.add_parser(
        "taxon-ids", help="Read and validate taxon IDs"
    )
    taxon_parser.add_argument("--file", required=True, help="Path to taxon ID file")

    # protein comp
    comp_parser = protein_sub.add_parser(
        "comp", help="Amino acid composition from FASTA"
    )
    comp_parser.add_argument("--fasta", required=True, help="Path to FASTA file")

    # protein rmsd-ca
    rmsd_parser = protein_sub.add_parser(
        "rmsd-ca", help="RMSD between CA atoms of two PDB files"
    )
    rmsd_parser.add_argument("--pdb-a", required=True, help="Path to first PDB file")
    rmsd_parser.add_argument("--pdb-b", required=True, help="Path to second PDB file")

    # Quality subcommands
    quality_parser = subparsers.add_parser("quality", help="Quality control commands")
    quality_sub = quality_parser.add_subparsers(dest="quality_command")

    batch_parser = quality_sub.add_parser(
        "batch-detect", help="Detect batch effects in a dataset"
    )
    batch_parser.add_argument(
        "--data", required=True, help="Path to CSV data matrix (samples × features)"
    )
    batch_parser.add_argument(
        "--batches", required=True, help="Path to batch labels file (one per line)"
    )
    batch_parser.add_argument(
        "--alpha", type=float, default=0.05, help="Significance threshold"
    )

    quality_run = quality_sub.add_parser(
        "run",
        help="Run the quality workflow: cross-code verification of docs against source",
    )
    quality_run.add_argument(
        "--output",
        type=Path,
        default=Path("output") / "cross_code_verification_report.md",
        help="Path for the Markdown verification report",
    )
    quality_run.add_argument(
        "--docs-dir",
        type=Path,
        default=Path("docs"),
        help="Documentation directory to verify",
    )
    quality_run.add_argument(
        "--src-dir", type=Path, default=Path("src"), help="Source directory to index"
    )
    quality_run.add_argument(
        "--include-historical",
        action="store_true",
        help="Include historical audit and validation snapshots",
    )
    quality_run.add_argument(
        "--strict-optional-imports",
        action="store_true",
        help="Treat optional third-party imports as violations",
    )
    quality_run.add_argument(
        "--verbose", action="store_true", help="Verbose verification logging"
    )

    # RNA subcommands
    rna_parser = subparsers.add_parser("rna", help="RNA-seq analysis commands")
    rna_sub = rna_parser.add_subparsers(dest="rna_command")

    rna_sub.add_parser("info", help="Show module capabilities")

    # GWAS subcommands
    gwas_parser = subparsers.add_parser("gwas", help="GWAS analysis commands")
    gwas_sub = gwas_parser.add_subparsers(dest="gwas_command")

    gwas_sub.add_parser("info", help="Show module capabilities")

    # gwas run subcommand
    gwas_run_parser = gwas_sub.add_parser("run", help="Run complete GWAS workflow")
    gwas_run_parser.add_argument(
        "--config", required=True, help="Path to GWAS configuration file (YAML/JSON)"
    )
    gwas_run_parser.add_argument(
        "--check", action="store_true", help="Validate configuration without executing"
    )
    gwas_run_parser.add_argument("--output-dir", help="Override output directory")

    # Life events subcommands
    life_parser = subparsers.add_parser(
        "life-events", help="Life event workflow commands"
    )
    life_sub = life_parser.add_subparsers(dest="life_events_command")

    life_predict = life_sub.add_parser(
        "predict", help="Predict outcomes for life event sequences"
    )
    life_predict.add_argument(
        "--events", required=True, help="Path to event sequences JSON"
    )
    life_predict.add_argument(
        "--model", required=True, help="Path to trained life-events model"
    )
    life_predict.add_argument(
        "--output", required=True, help="Output directory for predictions"
    )

    life_interpret = life_sub.add_parser(
        "interpret", help="Create a life-events interpretation report"
    )
    life_interpret.add_argument(
        "--model", required=True, help="Path to trained life-events model"
    )
    life_interpret.add_argument(
        "--sequences", required=True, help="Path to event sequences JSON"
    )
    life_interpret.add_argument(
        "--output", required=True, help="Output directory for report"
    )

    # Simulation subcommands
    simulation_parser = subparsers.add_parser(
        "simulation", help="Synthetic data simulation commands"
    )
    simulation_sub = simulation_parser.add_subparsers(dest="simulation_command")
    simulation_run = simulation_sub.add_parser(
        "run", help="Run a simulation workflow and save its result JSON"
    )
    simulation_run.add_argument(
        "--model",
        default="sequence_evolution",
        help=(
            "Simulation type: sequence_evolution, population_genetics, rna_expression, "
            "agent_ecosystem, predator_prey or competition"
        ),
    )
    simulation_run.add_argument(
        "--n",
        type=int,
        default=None,
        help="Simulation size override (population size, agent count or sample count depending on --model)",
    )
    simulation_run.add_argument(
        "--output",
        default="output/simulation",
        help="Output directory for the simulation result JSON",
    )

    # Ontology subcommands
    ontology_parser = subparsers.add_parser(
        "ontology", help="Ontology analysis workflow commands"
    )
    ontology_sub = ontology_parser.add_subparsers(dest="ontology_command")
    ontology_run = ontology_sub.add_parser(
        "run", help="Run the GO/HPO ontology enrichment workflow (stage 10)"
    )
    ontology_run.add_argument(
        "--input", required=True, help="Path to the workflow YAML configuration"
    )
    ontology_run.add_argument(
        "--phenotype",
        required=True,
        help="Phenotype label for the results subdirectory",
    )
    ontology_run.add_argument(
        "--model", required=True, help="Model label for the results subdirectory"
    )

    # Phenotype subcommands
    phenotype_parser = subparsers.add_parser(
        "phenotype", help="Phenotype analysis pipeline commands"
    )
    phenotype_sub = phenotype_parser.add_subparsers(dest="phenotype_command")
    phenotype_run = phenotype_sub.add_parser(
        "run", help="Run a phenotype analysis pipeline over a JSON dataset"
    )
    phenotype_run.add_argument(
        "--input", required=True, help="Path to phenotype data JSON (list of records)"
    )
    phenotype_run.add_argument(
        "--type",
        default="morphological",
        choices=("morphological", "behavioral", "chemical", "electronic", "sonic"),
        help="Phenotype domain to analyze",
    )
    phenotype_run.add_argument(
        "--output",
        default="output/phenotype",
        help="Output directory for the pipeline result JSON",
    )

    # Networks subcommands
    networks_parser = subparsers.add_parser(
        "networks", help="Network analysis workflow commands"
    )
    networks_sub = networks_parser.add_subparsers(dest="networks_command")
    networks_run = networks_sub.add_parser(
        "run", help="Build a network from an edge list and analyze it"
    )
    networks_run.add_argument(
        "--input",
        required=True,
        help="Path to an edge list CSV with 'source' and 'target' columns (optional 'weight')",
    )
    networks_run.add_argument(
        "--output",
        default="output/networks",
        help="Output directory for the exported network and metrics",
    )

    args = parser.parse_args()

    if args.modules:
        _list_modules()
        return 0

    if args.command == "protein":
        return _handle_protein(args)

    if args.command == "quality":
        return _handle_quality(args)

    if args.command == "rna":
        return _handle_rna(args)

    if args.command == "gwas":
        return _handle_gwas(args)

    if args.command == "life-events":
        return _handle_life_events(args)

    if args.command == "simulation":
        return _handle_simulation(args)

    if args.command == "ontology":
        return _handle_ontology(args)

    if args.command == "phenotype":
        return _handle_phenotype(args)

    if args.command == "networks":
        return _handle_networks(args)

    # If no arguments provided, show help
    if len(sys.argv) == 1:
        parser.print_help()
        return 0

    # A parsed command was not handled; never exit 0 without output.
    parser.print_help(sys.stderr)
    return 1


# Mapping from simulation type to the SimulationConfig attribute that --n
# overrides. Each workflow has one canonical "size" parameter.
_SIMULATION_SIZE_ATTRS = {
    "sequence_evolution": "population_size",
    "population_genetics": "population_size",
    "rna_expression": "n_samples",
    "agent_ecosystem": "n_agents",
    "predator_prey": "n_agents",
    "competition": "n_agents",
}


def _handle_protein(args: argparse.Namespace) -> int:
    """Handle protein subcommands."""
    import numpy as np

    cmd = args.protein_command

    if cmd == "taxon-ids":
        from .protein.sequence.proteomes import read_taxon_ids

        ids = read_taxon_ids(Path(args.file))
        print(" ".join(str(taxon_id) for taxon_id in ids))
        return 0

    elif cmd == "comp":
        from .protein.sequence.sequences import amino_acid_composition, read_fasta

        sequences = read_fasta(Path(args.fasta))
        for name, seq in sequences.items():
            comp = amino_acid_composition(seq)
            parts = [
                f"{aa}:{frac:.4f}" for aa, frac in sorted(comp.items()) if frac > 0
            ]
            print(f"{name}\t{','.join(parts)}")
        return 0

    elif cmd == "rmsd-ca":
        from .protein.structure.general import compute_rmsd_kabsch
        from .protein.structure.io import read_pdb_ca_coordinates

        ca_a = read_pdb_ca_coordinates(Path(args.pdb_a))
        ca_b = read_pdb_ca_coordinates(Path(args.pdb_b))
        rmsd = compute_rmsd_kabsch(np.array(ca_a), np.array(ca_b))
        print(f"{rmsd:.6f}")
        return 0

    print("Error: unknown or missing protein subcommand. See --help.", file=sys.stderr)
    return 1


def _handle_quality(args: argparse.Namespace) -> int:
    """Handle quality subcommands."""
    import numpy as np

    cmd = args.quality_command

    if cmd == "batch-detect":
        from .quality.batch.detection import detect_batch_effects

        data = np.loadtxt(args.data, delimiter=",", skiprows=1)
        batch_labels = Path(args.batches).read_text().strip().split("\n")
        report = detect_batch_effects(data, batch_labels, alpha=args.alpha)
        print(f"Samples: {report.n_samples}, Batches: {report.n_batches}")
        print(f"Batch variance: {report.pvca_variance['batch']:.3f}")
        print(f"Silhouette score: {report.silhouette_score:.3f}")
        print(f"Severity: {report.severity}")
        print(f"Significant features: {report.n_significant_features}")
        return 0

    if cmd == "run":
        from metainformant.quality.doc_verification import run as run_quality_workflow

        workflow_args = argparse.Namespace(
            verbose=args.verbose,
            docs_dir=args.docs_dir,
            src_dir=args.src_dir,
            output=args.output,
            include_historical=args.include_historical,
            strict_optional_imports=args.strict_optional_imports,
        )
        return run_quality_workflow(workflow_args)

    print("Error: unknown or missing quality subcommand. See --help.", file=sys.stderr)
    return 1


def _handle_rna(args: argparse.Namespace) -> int:
    """Handle RNA subcommands."""
    cmd = args.rna_command

    if cmd == "info":
        print("RNA-seq Analysis Module")
        print(
            "Sub-packages: amalgkit, analysis, core, deconvolution, engine, retrieval, splicing"
        )
        print("Import: from metainformant import rna")
        return 0

    print("Error: unknown or missing rna subcommand. See --help.", file=sys.stderr)
    return 1


def _handle_gwas(args: argparse.Namespace) -> int:
    """Handle GWAS subcommands."""
    cmd = args.gwas_command

    if cmd == "info":
        print("GWAS Analysis Module")
        print("Sub-packages: analysis, data, finemapping, heritability, visualization")
        print("Import: from metainformant import gwas")
        print("\nCLI Command: python -m metainformant gwas run --config <config.yaml>")
        return 0

    elif cmd == "run":
        # Import here to avoid heavy dependencies unless used
        try:
            from metainformant.gwas.workflow.workflow_execution import (
                execute_gwas_workflow,
            )
        except ImportError as e:
            print(f"Error: GWAS module dependencies not available: {e}")
            return 1

        config_path = args.config
        check_mode = args.check or args.config is None  # If no config, treat as check

        if not Path(config_path).exists():
            print(f"Error: Configuration file not found: {config_path}")
            return 1

        try:
            # Load configuration
            from metainformant.gwas.workflow.workflow_config import load_gwas_config

            config = load_gwas_config(config_path)

            # Override output directory if specified
            if args.output_dir:
                config["output_dir"] = args.output_dir

            # Execute or check
            if check_mode:
                result = execute_gwas_workflow(config, check=True)
                if result.get("status") == "validated":
                    print("✓ Configuration is valid")
                    return 0
                else:
                    print("✗ Configuration validation failed:")
                    for err in result.get("errors", []):
                        print(f"  - {err}")
                    return 1
            else:
                print(f"Starting GWAS workflow with config: {config_path}")
                result = execute_gwas_workflow(config, check=False)

                if result.get("success"):
                    print("✓ GWAS workflow completed successfully")
                    output_dir = result.get("output_dir", ".")
                    print(f"Results saved to: {output_dir}")
                    return 0
                else:
                    print("✗ GWAS workflow failed")
                    for err in result.get("errors", []):
                        print(f"  - {err}")
                    return 1

        except Exception as e:
            print(f"Error executing GWAS workflow: {e}")
            import traceback

            traceback.print_exc()
            return 1

    print("Error: unknown or missing gwas subcommand. See --help.", file=sys.stderr)
    return 1


def _load_life_event_sequences(path: str | Path) -> list:
    """Load EventSequence records from JSON for CLI helpers."""
    from metainformant.core.io.io import load_json
    from metainformant.life_events.core.events import EventSequence

    data = load_json(path)
    if isinstance(data, dict) and "sequences" in data:
        data = data["sequences"]
    if isinstance(data, dict):
        data = [data]
    return [EventSequence.from_dict(item) for item in data]


def _handle_life_events(args: argparse.Namespace) -> int:
    """Handle life-events subcommands."""
    from metainformant.core.io.io import dump_json
    from metainformant.life_events.core.utils import convert_sequences_to_tokens
    from metainformant.life_events.models.predictor import EventSequencePredictor

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.life_events_command == "predict":
        sequences = _load_life_event_sequences(args.events)
        predictor = EventSequencePredictor.load_model(args.model)
        tokens = convert_sequences_to_tokens(sequences)
        predictions = predictor.predict(tokens)
        prediction_values = (
            predictions.tolist()
            if hasattr(predictions, "tolist")
            else list(predictions)
        )
        probabilities = None
        if predictor.task_type == "classification":
            try:
                probabilities = predictor.predict_proba(tokens)
            except (AttributeError, ValueError):
                probabilities = None

        entries = []
        for i, (sequence, prediction) in enumerate(zip(sequences, prediction_values)):
            entry = {"sequence_id": sequence.person_id, "prediction": prediction}
            if probabilities is not None:
                raw_classes = getattr(predictor, "classes_", None)
                classes = (
                    raw_classes.tolist()
                    if raw_classes is not None and hasattr(raw_classes, "tolist")
                    else raw_classes
                )
                assert classes is not None
                prob_row = probabilities[i]
                if prob_row is None or not hasattr(prob_row, "tolist"):
                    continue
                entry["probabilities"] = {
                    str(cls): float(prob)
                    for cls, prob in zip(classes, prob_row.tolist())
                }
            entries.append(entry)

        payload = {
            "n_sequences": len(sequences),
            "model_path": str(args.model),
            "model_type": predictor.model_type,
            "task_type": predictor.task_type,
            "predictions": entries,
        }
        if predictor.task_type == "regression" and prediction_values:
            import numpy as np

            values = np.asarray(prediction_values, dtype=float)
            statistics = {
                "mean": float(values.mean()),
                "min": float(values.min()),
                "max": float(values.max()),
            }
            payload["statistics"] = statistics
            print(f"Mean: {statistics['mean']:.6f}")
        dump_json(payload, output_dir / "predictions.json")
        return 0

    if args.life_events_command == "interpret":
        sequences = _load_life_event_sequences(args.sequences)
        predictor = EventSequencePredictor.load_model(args.model)
        tokens = convert_sequences_to_tokens(sequences)
        predictions = predictor.predict(tokens)
        from metainformant.life_events.analysis.interpretability import (
            event_importance,
            feature_attribution,
            temporal_patterns,
        )

        embeddings = predictor.embeddings
        try:
            importance = event_importance(
                predictor, tokens, embeddings, method="permutation"
            )
        except ValueError:
            importance = event_importance(tokens)
        temporal = temporal_patterns(tokens, predictions)
        try:
            attribution = feature_attribution(predictor, tokens, embeddings)
        except ValueError:
            attribution = {"attributions": {}, "method": "unavailable"}

        report = {
            "model_path": str(args.model),
            "n_sequences": len(sequences),
            "model_type": predictor.model_type,
            "task_type": predictor.task_type,
            "predictions": predictions.tolist()
            if hasattr(predictions, "tolist")
            else list(predictions),
            "interpretations": {
                "event_importance": importance,
                "temporal_patterns": temporal,
                "feature_attribution": attribution,
            },
        }
        dump_json(report, output_dir / "interpretation_report.json")
        return 0

    print(
        "Error: unknown or missing life-events subcommand. See --help.", file=sys.stderr
    )
    return 1


def _handle_simulation(args: argparse.Namespace) -> int:
    """Handle simulation subcommands."""
    if args.simulation_command != "run":
        print(
            "Error: unknown or missing simulation subcommand. See --help.",
            file=sys.stderr,
        )
        return 1

    from metainformant.core.utils.errors import ConfigError, ValidationError
    from metainformant.simulation.workflow.workflow import (
        SimulationConfig,
        run_simulation_workflow,
    )

    config_kwargs: dict[str, Any] = {
        "simulation_type": args.model,
        "output_dir": args.output,
    }
    size_attr = _SIMULATION_SIZE_ATTRS.get(args.model)
    if args.n is not None and size_attr is not None:
        config_kwargs[size_attr] = args.n

    try:
        config = SimulationConfig(**config_kwargs)
    except (ValidationError, ConfigError, ValueError) as exc:
        print(f"Error: invalid simulation configuration: {exc}", file=sys.stderr)
        return 1

    result = run_simulation_workflow(config)
    print(f"Simulation type: {config.simulation_type}")
    output_file = result.get("output_file")
    if output_file:
        print(f"Result saved to: {output_file}")
    return 0


def _handle_ontology(args: argparse.Namespace) -> int:
    """Handle ontology subcommands."""
    if args.ontology_command != "run":
        print(
            "Error: unknown or missing ontology subcommand. See --help.",
            file=sys.stderr,
        )
        return 1

    config_path = Path(args.input)
    if not config_path.exists():
        print(f"Error: workflow config not found: {config_path}", file=sys.stderr)
        return 1

    from metainformant.ontology.workflow.run_ontology import run_ontology_analysis

    try:
        import yaml

        with open(config_path) as f:
            config = yaml.safe_load(f) or {}
    except Exception as exc:
        print(f"Error: could not load workflow config: {exc}", file=sys.stderr)
        return 1
    if not isinstance(config, dict):
        print("Error: workflow config must be a YAML mapping", file=sys.stderr)
        return 1

    # run_ontology_analysis exits with its own status code when stage inputs
    # are missing, so let SystemExit propagate from main().
    print(f"Ontology workflow: phenotype={args.phenotype} model={args.model}")
    run_ontology_analysis(config, args.phenotype, args.model, Path.cwd())
    return 0


def _handle_phenotype(args: argparse.Namespace) -> int:
    """Handle phenotype subcommands."""
    if args.phenotype_command != "run":
        print(
            "Error: unknown or missing phenotype subcommand. See --help.",
            file=sys.stderr,
        )
        return 1

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: phenotype data not found: {input_path}", file=sys.stderr)
        return 1

    from metainformant.phenotype.workflow.pipeline import (
        PhenotypePipeline,
        PipelineConfig,
    )

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    config = PipelineConfig(phenotype_types=[args.type], input_path=str(input_path))
    result = PhenotypePipeline(config).run()

    result_path = output_dir / "pipeline_result.json"
    result.save_json(result_path)
    print(f"Pipeline: {config.name} (type: {args.type})")
    for step_name, step_output in result.outputs.items():
        status = (
            step_output.get("status", "done")
            if isinstance(step_output, dict)
            else "done"
        )
        print(f"  {step_name}: {status}")
    print(f"Result saved to: {result_path}")
    for error in result.errors:
        print(f"Error: {error}", file=sys.stderr)
    return 0 if result.success else 1


def _handle_networks(args: argparse.Namespace) -> int:
    """Handle networks subcommands."""
    if args.networks_command != "run":
        print(
            "Error: unknown or missing networks subcommand. See --help.",
            file=sys.stderr,
        )
        return 1

    import csv

    from metainformant.networks.config.config import NetworkWorkflowConfig
    from metainformant.networks.workflow.workflow import NetworkWorkflow

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: edge list not found: {input_path}", file=sys.stderr)
        return 1

    edges: list[tuple[str, str] | tuple[str, str, float]] = []
    with open(input_path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        if "source" not in fieldnames or "target" not in fieldnames:
            print(
                "Error: edge list CSV must contain 'source' and 'target' columns",
                file=sys.stderr,
            )
            return 1
        for row in reader:
            source = (row.get("source") or "").strip()
            target = (row.get("target") or "").strip()
            if not source or not target:
                continue
            weight_raw = (row.get("weight") or "").strip()
            if weight_raw:
                try:
                    edges.append((source, target, float(weight_raw)))
                except ValueError:
                    print(
                        f"Error: invalid weight for edge {source}->{target}: {weight_raw!r}",
                        file=sys.stderr,
                    )
                    return 1
            else:
                edges.append((source, target))

    if not edges:
        print("Error: no edges found in edge list", file=sys.stderr)
        return 1

    output_dir = Path(args.output)
    workflow = NetworkWorkflow(NetworkWorkflowConfig(output_dir=str(output_dir)))
    workflow.build_network(edges=edges).detect_communities().analyze_metrics()
    exported = workflow.export_results(str(output_dir))
    summary = workflow.summary()
    network_summary = summary.get("network", {})
    print(
        f"Network: {network_summary.get('n_nodes', 0)} nodes, {network_summary.get('n_edges', 0)} edges"
    )
    print(f"Communities: {(summary.get('communities') or {}).get('n_communities', 0)}")
    print(f"Results exported to: {output_dir} ({len(exported)} files)")
    return 0


def _list_modules() -> None:
    """List all available modules."""
    modules = [
        ("core", "Shared utilities and infrastructure"),
        ("dna", "DNA sequence analysis and genomics"),
        ("rna", "RNA-seq workflows and amalgkit integration"),
        ("protein", "Protein sequence and structure analysis"),
        ("gwas", "Genome-wide association studies"),
        ("math", "Mathematical biology and theoretical modeling"),
        ("information", "Information-theoretic analysis"),
        ("life_events", "Life course and temporal analysis"),
        ("visualization", "Plotting and visualization tools"),
        ("networks", "Biological network analysis"),
        ("multiomics", "Cross-omics data integration"),
        ("singlecell", "Single-cell RNA-seq analysis"),
        ("simulation", "Synthetic data generation"),
        ("quality", "Data quality control"),
        ("ml", "Machine learning for biological data"),
        ("ontology", "Gene ontology and functional annotation"),
        ("phenotype", "Phenotypic trait analysis"),
        ("ecology", "Ecological and community analysis"),
        ("epigenome", "Epigenomic data analysis"),
        ("longread", "Long-read sequencing analysis (ONT, PacBio)"),
        ("structural_variants", "Structural variant detection and analysis"),
        ("spatial", "Spatial transcriptomics analysis"),
        ("metagenomics", "Microbiome and metagenomic analysis"),
        ("pharmacogenomics", "Clinical pharmacogenomic variant analysis"),
        ("metabolomics", "Metabolite identification and pathway analysis"),
        ("cloud", "Cloud deployment helpers and GCP workflow utilities"),
        (
            "mcp",
            "MCP stdio JSON-RPC 2.0 server (`python -m metainformant.mcp.server`) and tool registry",
        ),
        ("menu", "Interactive menu and discovery system"),
    ]

    print("Available METAINFORMANT modules:")
    print("=" * 50)

    for name, description in modules:
        print(f"  {name:15} - {description}")

    print("\nImport modules in Python:")
    print("  from metainformant import dna, rna, protein  # etc.")


if __name__ == "__main__":
    sys.exit(main())
