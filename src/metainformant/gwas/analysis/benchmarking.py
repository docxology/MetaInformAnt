"""GWAS compute-time benchmarking and runtime extrapolation.

This module provides tools for:
- Timing pilot GWAS runs on data subsets (limited samples or loci)
- Extrapolating full-genome compute times from pilot timings using
  known computational complexity models for each pipeline step
- Generating human-readable runtime estimates

Scaling models used:
    QC / MAF filtering:   O(n · m)   — linear in samples × variants
    LD pruning:           O(w² · m)  — quadratic in window, linear in variants
    PCA:                  O(m · k²)  — linear in variants, quadratic in components
    Kinship (GRM):        O(n² · m)  — quadratic in samples, linear in variants
    Association testing:  O(n · m)   — linear in samples × variants
    Fine-mapping (SuSiE): O(k³)      — cubic in credible-set region size
    Heritability (LDSC):  O(m)       — linear in variants
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class StepTiming:
    """Timing result for a single GWAS pipeline step."""

    step_name: str
    elapsed_seconds: float
    n_samples: int
    n_variants: int
    extra_params: Dict[str, Any] = field(default_factory=dict)

    @property
    def throughput(self) -> float:
        """Tests per second (n_samples × n_variants / elapsed)."""
        total = self.n_samples * self.n_variants
        return total / self.elapsed_seconds if self.elapsed_seconds > 0 else 0.0


@dataclass
class ComputeTimeEstimate:
    """Full-genome runtime estimate extrapolated from pilot timings."""

    total_seconds: float
    total_human: str
    per_step: Dict[str, float]
    pilot_n_samples: int
    pilot_n_variants: int
    target_n_samples: int
    target_n_variants: int
    scaling_factors: Dict[str, float] = field(default_factory=dict)

    def summary(self) -> str:
        """Return a human-readable summary."""
        lines = [
            f"Estimated total runtime: {self.total_human}",
            f"  Pilot: {self.pilot_n_samples} samples × {self.pilot_n_variants} variants",
            f"  Target: {self.target_n_samples} samples × {self.target_n_variants} variants",
            "",
            "Per-step estimates:",
        ]
        for step, secs in sorted(self.per_step.items()):
            hrs = secs / 3600
            factor = self.scaling_factors.get(step, 1.0)
            lines.append(f"  {step:30s}  {_format_duration(secs):>12s}  (×{factor:.1f})")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Scaling models
# ---------------------------------------------------------------------------

# Each model maps (pilot_n, pilot_m, target_n, target_m) → scale factor
# where n = samples, m = variants.

SCALING_MODELS: Dict[str, str] = {
    "parse_vcf":              "n_m",       # O(n·m)
    "qc_filters":             "n_m",       # O(n·m)
    "maf_filter":             "n_m",       # O(n·m)
    "hwe_test":               "n_m",       # O(n·m)
    "ld_pruning":             "m",         # O(w²·m) — window fixed, linear in m
    "population_structure":   "n2_k2",     # O(n²·m) dominated by kinship
    "pca":                    "m_k2",      # O(m·k²)
    "kinship":                "n2_m",      # O(n²·m)
    "association_testing":    "n_m",       # O(n·m)
    "multiple_testing":       "m",         # O(m)
    "fine_mapping":           "k3",        # O(k³) per region
    "heritability":           "m",         # O(m) for LDSC
    "visualization":          "m",         # O(m) for Manhattan/QQ
    "summary_stats":          "m",         # O(m)
    "annotation":             "m",         # O(m)
}


def scaling_model(
    model: str,
    pilot_n: int,
    pilot_m: int,
    target_n: int,
    target_m: int,
    *,
    k_pilot: int = 10,
    k_target: int = 10,
) -> float:
    """Compute the scaling factor for a given complexity model.

    Args:
        model: Complexity model identifier (e.g. "n_m", "n2_m", "m").
        pilot_n: Number of samples in pilot run.
        pilot_m: Number of variants in pilot run.
        target_n: Number of samples in target (full) run.
        target_m: Number of variants in target (full) run.
        k_pilot: Number of PCA components / fine-mapping region size in pilot.
        k_target: Number of PCA components / fine-mapping region size in target.

    Returns:
        Multiplicative scaling factor (target_time ≈ pilot_time × factor).

    Raises:
        ValueError: If pilot dimensions are zero.
    """
    if pilot_n <= 0 or pilot_m <= 0:
        raise ValueError(f"Pilot dimensions must be positive: n={pilot_n}, m={pilot_m}")

    if model == "n_m":
        # O(n · m)
        return (target_n * target_m) / (pilot_n * pilot_m)

    elif model == "n2_m":
        # O(n² · m)
        return (target_n ** 2 * target_m) / (pilot_n ** 2 * pilot_m)

    elif model == "m":
        # O(m)
        return target_m / pilot_m

    elif model == "m_k2":
        # O(m · k²)
        return (target_m * k_target ** 2) / (pilot_m * k_pilot ** 2)

    elif model == "n2_k2":
        # O(n² · m) dominated term
        return (target_n ** 2 * target_m) / (pilot_n ** 2 * pilot_m)

    elif model == "k3":
        # O(k³) — region-level, approximate with variant ratio
        return (target_m / pilot_m) ** 1.5  # sub-cubic approximation

    else:
        logger.warning(f"Unknown scaling model '{model}', using linear (O(m)) fallback")
        return target_m / pilot_m


# ---------------------------------------------------------------------------
# Pilot-run benchmarking
# ---------------------------------------------------------------------------


def benchmark_subset_run(
    vcf_path: Union[str, Path],
    phenotype_path: Union[str, Path],
    config: Dict[str, Any],
    *,
    max_samples: Optional[int] = None,
    max_variants: Optional[int] = None,
) -> List[StepTiming]:
    """Run a GWAS pipeline on a data subset and time each step.

    This executes the real ``run_gwas`` pipeline on a subset of data,
    capturing the wall-clock time for each pipeline step.

    Args:
        vcf_path: Path to VCF file.
        phenotype_path: Path to phenotype file.
        config: GWAS configuration dictionary.
        max_samples: Limit the number of samples (None = use all).
        max_variants: Limit the number of variants (None = use all).

    Returns:
        List of StepTiming objects, one per pipeline step.
    """
    from metainformant.gwas.analysis.quality import parse_vcf_full
    from metainformant.gwas.workflow.workflow_execution import run_gwas

    vcf_path = Path(vcf_path)
    phenotype_path = Path(phenotype_path)

    # --- Determine actual data dimensions ---
    t0 = time.perf_counter()
    vcf_data = parse_vcf_full(vcf_path)
    parse_time = time.perf_counter() - t0

    actual_n_samples = len(vcf_data.get("samples", []))
    actual_n_variants = len(vcf_data.get("variants", []))

    n_samples_used = min(max_samples, actual_n_samples) if max_samples else actual_n_samples
    n_variants_used = min(max_variants, actual_n_variants) if max_variants else actual_n_variants

    logger.info(
        f"Benchmarking subset: {n_samples_used}/{actual_n_samples} samples, "
        f"{n_variants_used}/{actual_n_variants} variants"
    )

    # --- Run the full pipeline and capture overall timing ---
    import tempfile

    with tempfile.TemporaryDirectory(prefix="gwas_bench_") as tmpdir:
        pilot_config = dict(config)
        if max_samples:
            pilot_config["max_samples"] = max_samples
        if max_variants:
            pilot_config["max_variants"] = max_variants

        t_start = time.perf_counter()
        result = run_gwas(vcf_path, phenotype_path, pilot_config, output_dir=tmpdir)
        total_time = time.perf_counter() - t_start

    # --- Extract per-step timings from result if available ---
    timings: List[StepTiming] = []
    steps_completed = result.get("steps_completed", [])

    # If the pipeline reports step_timings, use them
    step_timings_raw = result.get("step_timings", {})

    if step_timings_raw:
        for step_name, elapsed in step_timings_raw.items():
            timings.append(StepTiming(
                step_name=step_name,
                elapsed_seconds=elapsed,
                n_samples=n_samples_used,
                n_variants=n_variants_used,
            ))
    else:
        # Fallback: report total time split proportionally by known complexity ratios
        proportions = _estimate_step_proportions(steps_completed)
        for step_name, proportion in proportions.items():
            timings.append(StepTiming(
                step_name=step_name,
                elapsed_seconds=total_time * proportion,
                n_samples=n_samples_used,
                n_variants=n_variants_used,
            ))

    # Add the VCF parse timing we measured directly
    if not any(t.step_name == "parse_vcf" for t in timings):
        timings.insert(0, StepTiming(
            step_name="parse_vcf",
            elapsed_seconds=parse_time,
            n_samples=n_samples_used,
            n_variants=n_variants_used,
        ))

    logger.info(
        f"Benchmark complete: {len(timings)} steps timed, total {total_time:.2f}s"
    )
    return timings


def extrapolate_full_genome_time(
    pilot_timings: List[StepTiming],
    target_n_samples: int,
    target_n_variants: int,
    *,
    custom_models: Optional[Dict[str, str]] = None,
) -> ComputeTimeEstimate:
    """Extrapolate full-genome runtime from pilot-run timings.

    Uses known computational complexity models to scale pilot-run timings
    to the target data dimensions.

    Args:
        pilot_timings: StepTiming list from ``benchmark_subset_run()``.
        target_n_samples: Number of samples in full dataset.
        target_n_variants: Number of variants in full dataset.
        custom_models: Optional mapping of step_name → scaling model to override
            the defaults in ``SCALING_MODELS``.

    Returns:
        ComputeTimeEstimate with per-step and total estimates.

    Example:
        >>> timings = benchmark_subset_run("small.vcf", "pheno.tsv", config)
        >>> est = extrapolate_full_genome_time(timings, 5000, 10_000_000)
        >>> print(est.summary())
    """
    if not pilot_timings:
        return ComputeTimeEstimate(
            total_seconds=0.0,
            total_human="0s",
            per_step={},
            pilot_n_samples=0,
            pilot_n_variants=0,
            target_n_samples=target_n_samples,
            target_n_variants=target_n_variants,
        )

    models = dict(SCALING_MODELS)
    if custom_models:
        models.update(custom_models)

    pilot_n = pilot_timings[0].n_samples
    pilot_m = pilot_timings[0].n_variants

    per_step: Dict[str, float] = {}
    factors: Dict[str, float] = {}

    for timing in pilot_timings:
        model_name = models.get(timing.step_name, "n_m")
        factor = scaling_model(
            model_name,
            pilot_n=timing.n_samples,
            pilot_m=timing.n_variants,
            target_n=target_n_samples,
            target_m=target_n_variants,
        )
        estimated_seconds = timing.elapsed_seconds * factor
        per_step[timing.step_name] = estimated_seconds
        factors[timing.step_name] = factor

    total = sum(per_step.values())

    return ComputeTimeEstimate(
        total_seconds=total,
        total_human=_format_duration(total),
        per_step=per_step,
        pilot_n_samples=pilot_n,
        pilot_n_variants=pilot_m,
        target_n_samples=target_n_samples,
        target_n_variants=target_n_variants,
        scaling_factors=factors,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _format_duration(seconds: float) -> str:
    """Format seconds into human-readable duration string.

    Args:
        seconds: Duration in seconds.

    Returns:
        Human-readable string like "2h 15m 30s" or "45.2s".
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        m, s = divmod(seconds, 60)
        return f"{int(m)}m {int(s)}s"
    else:
        h, remainder = divmod(seconds, 3600)
        m, s = divmod(remainder, 60)
        return f"{int(h)}h {int(m)}m {int(s)}s"


def _estimate_step_proportions(steps: List[str]) -> Dict[str, float]:
    """Estimate relative time proportions for steps when per-step timings unavailable.

    Based on empirical observations from real GWAS runs (Apis mellifera cohort).

    Args:
        steps: List of completed step names.

    Returns:
        Dictionary mapping step name to proportion of total time.
    """
    # Empirical proportions from Apis mellifera 188-sample, 3000-variant run
    default_proportions = {
        "parse_vcf": 0.05,
        "qc_filters": 0.05,
        "ld_pruning": 0.05,
        "population_structure": 0.15,
        "load_phenotypes": 0.01,
        "association_testing": 0.50,
        "multiple_testing_correction": 0.02,
        "fine_mapping": 0.07,
        "summary_stats": 0.02,
        "annotation": 0.03,
        "visualization": 0.05,
    }

    # Filter to steps that were actually completed
    matched = {s: default_proportions.get(s, 0.05) for s in steps}

    # Normalize to sum to 1.0
    total = sum(matched.values())
    if total > 0:
        matched = {k: v / total for k, v in matched.items()}

    return matched
