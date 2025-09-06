from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from metainformant.core import io as core_io
from metainformant.core import paths as core_paths

from .model import noise, simulate_generation, simulate_generations
from .plotting import display_generation, display_s_vs_q, display_sq_vs_w


def run_replay(dest_dir: Path) -> None:
    # Create outputs subfolder within selection_experiments directory
    outputs_dir = dest_dir / "outputs"
    core_io.ensure_directory(outputs_dir)
    np.random.seed(428)

    print(f"Running natural selection experiments...")
    print(f"Outputs will be saved to: {outputs_dir}")
    print()

    print("1. Simulating structural trait vs quality relationship...")
    r1 = simulate_generation(np.random.normal(0, 1, 1_000_000), s_hat=0.8)
    display_s_vs_q(r1, dest=outputs_dir / "plot-s-vs-q.png")
    print("   ✓ Generated plot-s-vs-q.png")

    print("2. Simulating traits vs fitness correlation...")
    r2 = simulate_generation(np.random.normal(0, 1, 1_000_000), s_hat=0.65)
    display_sq_vs_w(r2, dest=outputs_dir / "plot-sq-vs-w.png")
    print("   ✓ Generated plot-sq-vs-w.png")

    print("3. Running rebound selection experiment...")
    gr1 = simulate_generations(
        generations=50,
        n=10_000,
        s_hat=0.6,
        delta_fn=noise(0, -1, 4.7),
        fitness_fn=noise(0.1, 0, 10),
    )
    display_generation(gr1, dest=outputs_dir / "plot-ns-rebound.png")
    print("   ✓ Generated plot-ns-rebound.png")

    print("4. Running inverse selection experiment...")
    gr2 = simulate_generations(
        generations=50,
        n=1_000,
        s_hat=0.6,
        delta_fn=noise(0, 6, 10),
        fitness_fn=noise(-0.02, 0, 3),
    )
    display_generation(gr2, dest=outputs_dir / "plot-ns-inverse.png")
    print("   ✓ Generated plot-ns-inverse.png")

    print("5. Running neutral selection experiment...")
    gr3 = simulate_generations(
        generations=30,
        n=1_000,
        s_hat=0.1,
        delta_fn=noise(0, 0, 0.1),
        fitness_fn=noise(1.6, 0, 1.2),
    )
    display_generation(gr3, dest=outputs_dir / "plot-ns.png")
    print("   ✓ Generated plot-ns.png")

    print("6. Running quality signal limit (QSL) experiment...")
    gr4 = simulate_generations(
        generations=20,
        n=1_000,
        s_hat=0.006,
        delta_fn=noise(0, 0, 0.04),
        fitness_fn=noise(0.6, 0, 0.1),
    )
    display_generation(gr4, dest=outputs_dir / "plot-ns-qsl.png")
    print("   ✓ Generated plot-ns-qsl.png")

    print()
    print(f"🎉 All experiments completed! Check {outputs_dir} for PNG outputs.")


def run_create_abstract(dest_dir: Path) -> None:
    """Create a composite graphical abstract by running all experiments and combining results."""
    from .plotting import create_composite_abstract

    print("Creating composite graphical abstract...")
    print(f"Output directory: {dest_dir}")
    print()

    # Run all experiments to gather results
    experiment_results = {}

    print("1. Running signal mapping experiment...")
    r1 = simulate_generation(np.random.normal(0, 1, 1_000_000), s_hat=0.8)
    experiment_results["signal_mapping"] = r1

    print("2. Running selection correlation experiment...")
    r2 = simulate_generation(np.random.normal(0, 1, 1_000_000), s_hat=0.65)
    experiment_results["selection_correlation"] = r2

    print("3. Running rebound evolution experiment...")
    gr1 = simulate_generations(
        generations=50,
        n=10_000,
        s_hat=0.6,
        delta_fn=noise(0, -1, 4.7),
        fitness_fn=noise(0.1, 0, 10),
    )
    experiment_results["rebound_evolution"] = gr1

    print("4. Running inverse evolution experiment...")
    gr2 = simulate_generations(
        generations=50,
        n=1_000,
        s_hat=0.6,
        delta_fn=noise(0, 6, 10),
        fitness_fn=noise(-0.02, 0, 3),
    )
    experiment_results["inverse_evolution"] = gr2

    print("5. Running neutral evolution experiment...")
    gr3 = simulate_generations(
        generations=30,
        n=1_000,
        s_hat=0.1,
        delta_fn=noise(0, 0, 0.1),
        fitness_fn=noise(1.6, 0, 1.2),
    )
    experiment_results["neutral_evolution"] = gr3

    print("6. Running quality signal limit experiment...")
    gr4 = simulate_generations(
        generations=20,
        n=1_000,
        s_hat=0.006,
        delta_fn=noise(0, 0, 0.04),
        fitness_fn=noise(0.6, 0, 0.1),
    )
    experiment_results["qsl_evolution"] = gr4

    print()
    print("Generating composite graphical abstract...")

    # Create the composite abstract
    abstract_path = dest_dir / "graphical_abstract.png"
    create_composite_abstract(experiment_results, dest=abstract_path)

    print("✅ Composite graphical abstract created!")
    print(f"📄 Saved to: {abstract_path}")
    print()
    print("The graphical abstract combines all experimental results into a single comprehensive figure,")
    print("showing the key findings from the natural selection experiments on signal processing.")


def add_math_selection_subparser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    sel_parser = subparsers.add_parser("selection", help="Natural selection experiments (from notebook)")
    sel_sub = sel_parser.add_subparsers(dest="sel_cmd")
    replay = sel_sub.add_parser("replay", help="Reproduce notebook figures and outputs")
    replay.add_argument(
        "--dest",
        default="src/metainformant/math/selection_experiments",
        help="Destination directory for outputs (outputs subfolder will be created)",
    )

    abstract = sel_sub.add_parser("abstract", help="Create composite graphical abstract")
    abstract.add_argument(
        "--dest", default="src/metainformant/math/selection_experiments", help="Destination directory for the abstract"
    )

    def _dispatch(args: argparse.Namespace) -> None:
        if args.sel_cmd == "replay":
            dest_dir = core_paths.expand_and_resolve(args.dest)
            run_replay(dest_dir)
        elif args.sel_cmd == "abstract":
            dest_dir = core_paths.expand_and_resolve(args.dest)
            run_create_abstract(dest_dir)

    sel_parser.set_defaults(func=_dispatch)
