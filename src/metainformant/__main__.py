from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Explicitly import the test runner module to avoid confusion with the repo's root-level tests/
from .tests.runner import run_all_tests
# defer DNA imports until the DNA subcommand is used to avoid optional deps at import time
from .rna.workflow import AmalgkitWorkflowConfig, execute_workflow
from .rna.configs import SpeciesProfile, AmalgkitRunLayout, build_step_params
from .rna.workflow import plan_workflow_with_params


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

    if args.command == "math":
        func = getattr(args, "func", None)
        if callable(func):
            func(args)
            return

    if args.command == "tests":
        exit_code = run_all_tests(args.pytest_args)
        sys.exit(exit_code)

    parser.print_help()


if __name__ == "__main__":
    main()


