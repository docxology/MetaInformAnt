"""Synthetic benchmark dataset generators for bioinformatics method evaluation.

Provides generators for classification, regression, clustering, and
differential expression benchmarks with known ground truth, as well as
synthetic GWAS and RNA-seq datasets. Includes evaluation functions and
a benchmark suite runner for comparing multiple methods.

All generators use pure Python with optional NumPy acceleration for
matrix operations.
"""

from __future__ import annotations

import math
import random
from typing import Any, Callable

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _set_seed(seed: int | None) -> None:
    """Set random seed for reproducibility."""
    if seed is not None:
        random.seed(seed)
        if HAS_NUMPY:
            np.random.seed(seed)


def _mean(values: list[float]) -> float:
    """Arithmetic mean."""
    return sum(values) / len(values) if values else 0.0


def _std(values: list[float]) -> float:
    """Sample standard deviation."""
    n = len(values)
    if n < 2:
        return 0.0
    mu = _mean(values)
    return math.sqrt(sum((x - mu) ** 2 for x in values) / (n - 1))


# ---------------------------------------------------------------------------
# General benchmark datasets
# ---------------------------------------------------------------------------


def generate_benchmark_dataset(
    task: str = "classification",
    n_samples: int = 100,
    n_features: int = 1000,
    difficulty: str = "medium",
    seed: int | None = None,
) -> dict:
    """Generate a synthetic benchmark dataset with known ground truth.

    Args:
        task: Benchmark task type: ``"classification"``, ``"regression"``,
            ``"clustering"``, or ``"de_genes"`` (differential expression).
        n_samples: Number of samples to generate.
        n_features: Number of features.
        difficulty: Difficulty level: ``"easy"``, ``"medium"``, or ``"hard"``.
            Controls signal-to-noise ratio and class separation.
        seed: Random seed for reproducibility.

    Returns:
        Dictionary with keys:
            - X: 2D list of shape ``(n_samples, n_features)``.
            - y: List of labels/values (class labels or continuous).
            - true_labels: Ground truth labels.
            - metadata: Dict with generation parameters.
            - task_description: Human-readable task description.

    Raises:
        ValueError: If task is unknown.
    """
    _set_seed(seed)

    noise_scale = {"easy": 0.3, "medium": 1.0, "hard": 3.0}.get(difficulty, 1.0)
    signal_scale = {"easy": 3.0, "medium": 1.5, "hard": 0.5}.get(difficulty, 1.5)

    if task == "classification":
        return _generate_classification(n_samples, n_features, noise_scale, signal_scale, seed)
    elif task == "regression":
        return _generate_regression(n_samples, n_features, noise_scale, signal_scale, seed)
    elif task == "clustering":
        return _generate_clustering(n_samples, n_features, noise_scale, signal_scale, seed)
    elif task == "de_genes":
        return _generate_de_genes(n_samples, n_features, noise_scale, signal_scale, seed)
    else:
        raise ValueError(f"Unknown task '{task}'. Use: classification, regression, clustering, de_genes")


def _generate_classification(
    n_samples: int,
    n_features: int,
    noise: float,
    signal: float,
    seed: int | None,
) -> dict:
    """Generate binary classification dataset."""
    n_informative = max(5, n_features // 20)
    X: list[list[float]] = []
    y: list[int] = []

    for i in range(n_samples):
        label = 0 if i < n_samples // 2 else 1
        row = [random.gauss(0, noise) for _ in range(n_features)]
        # Add signal to informative features
        for j in range(n_informative):
            row[j] += signal * (1 if label == 1 else -1)
        X.append(row)
        y.append(label)

    # Shuffle
    combined = list(zip(X, y))
    random.shuffle(combined)
    X = [c[0] for c in combined]
    y = [c[1] for c in combined]

    logger.info(
        "Generated classification dataset: %d samples, %d features, %d informative",
        n_samples,
        n_features,
        n_informative,
    )

    return {
        "X": X,
        "y": y,
        "true_labels": list(y),
        "metadata": {
            "task": "classification",
            "n_informative": n_informative,
            "noise_scale": noise,
            "signal_scale": signal,
            "seed": seed,
        },
        "task_description": (
            f"Binary classification with {n_informative} informative features " f"out of {n_features} total."
        ),
    }


def _generate_regression(
    n_samples: int,
    n_features: int,
    noise: float,
    signal: float,
    seed: int | None,
) -> dict:
    """Generate regression dataset."""
    n_informative = max(5, n_features // 20)
    true_coeffs = [random.gauss(0, signal) if j < n_informative else 0.0 for j in range(n_features)]

    X: list[list[float]] = []
    y: list[float] = []

    for _ in range(n_samples):
        row = [random.gauss(0, 1) for _ in range(n_features)]
        target = sum(row[j] * true_coeffs[j] for j in range(n_features))
        target += random.gauss(0, noise)
        X.append(row)
        y.append(target)

    return {
        "X": X,
        "y": y,
        "true_labels": y,
        "metadata": {
            "task": "regression",
            "n_informative": n_informative,
            "true_coefficients": true_coeffs[:n_informative],
            "noise_scale": noise,
            "seed": seed,
        },
        "task_description": (f"Regression with {n_informative} informative features."),
    }


def _generate_clustering(
    n_samples: int,
    n_features: int,
    noise: float,
    signal: float,
    seed: int | None,
) -> dict:
    """Generate clustering dataset with known clusters."""
    n_clusters = 3
    samples_per_cluster = n_samples // n_clusters

    X: list[list[float]] = []
    y: list[int] = []

    centers = [[signal * random.gauss(0, 1) for _ in range(min(10, n_features))] for _ in range(n_clusters)]

    for c in range(n_clusters):
        n_c = samples_per_cluster if c < n_clusters - 1 else n_samples - len(X)
        for _ in range(n_c):
            row = [random.gauss(0, noise) for _ in range(n_features)]
            for j in range(min(10, n_features)):
                row[j] += centers[c][j]
            X.append(row)
            y.append(c)

    combined = list(zip(X, y))
    random.shuffle(combined)
    X = [c[0] for c in combined]
    y = [c[1] for c in combined]

    return {
        "X": X,
        "y": y,
        "true_labels": list(y),
        "metadata": {
            "task": "clustering",
            "n_clusters": n_clusters,
            "noise_scale": noise,
            "signal_scale": signal,
            "seed": seed,
        },
        "task_description": f"Clustering with {n_clusters} true clusters.",
    }


def _generate_de_genes(
    n_samples: int,
    n_features: int,
    noise: float,
    signal: float,
    seed: int | None,
) -> dict:
    """Generate differential expression dataset with known DE genes."""
    n_de = max(10, n_features // 10)
    half = n_samples // 2

    X: list[list[float]] = []
    y: list[int] = []
    de_genes = list(range(n_de))
    fold_changes = [signal * (1 + random.random()) * random.choice([-1, 1]) for _ in range(n_de)]

    for i in range(n_samples):
        group = 0 if i < half else 1
        # Base expression ~ Poisson-like with log-normal noise
        row = [max(0.0, random.gauss(5.0, noise)) for _ in range(n_features)]
        if group == 1:
            for j, fc in zip(de_genes, fold_changes):
                row[j] = max(0.0, row[j] + fc)
        X.append(row)
        y.append(group)

    return {
        "X": X,
        "y": y,
        "true_labels": list(y),
        "metadata": {
            "task": "de_genes",
            "n_de_genes": n_de,
            "de_gene_indices": de_genes,
            "true_fold_changes": fold_changes,
            "seed": seed,
        },
        "task_description": (f"Differential expression with {n_de} true DE genes " f"out of {n_features} total."),
    }


# ---------------------------------------------------------------------------
# Synthetic GWAS data
# ---------------------------------------------------------------------------


def generate_synthetic_variants(
    n_variants: int = 1000,
    n_causal: int = 10,
    effect_sizes: list[float] | None = None,
    maf_range: tuple = (0.01, 0.5),
    seed: int | None = None,
    n_samples: int = 500,
) -> dict:
    """Generate synthetic GWAS data with known causal variants.

    Args:
        n_variants: Total number of variants.
        n_causal: Number of causal variants.
        effect_sizes: Optional list of effect sizes for causal variants.
            If ``None``, random effects are generated.
        maf_range: Tuple of (min_maf, max_maf) for allele frequencies.
        seed: Random seed.
        n_samples: Number of individuals.

    Returns:
        Dictionary with keys:
            - genotypes: 2D list ``(n_samples, n_variants)`` of 0/1/2 dosages.
            - phenotypes: List of continuous phenotype values.
            - causal_variants: List of causal variant indices.
            - true_effects: List of true effect sizes for causal variants.
            - heritability: Approximate heritability of the trait.
            - maf: List of minor allele frequencies.
    """
    _set_seed(seed)

    # Generate MAFs
    mafs = [random.uniform(maf_range[0], maf_range[1]) for _ in range(n_variants)]

    # Generate genotypes (Hardy-Weinberg)
    genotypes: list[list[int]] = []
    for _ in range(n_samples):
        row = []
        for maf in mafs:
            p = maf
            r = random.random()
            if r < (1 - p) ** 2:
                g = 0
            elif r < (1 - p) ** 2 + 2 * p * (1 - p):
                g = 1
            else:
                g = 2
            row.append(g)
        genotypes.append(row)

    # Select causal variants
    causal_indices = random.sample(range(n_variants), min(n_causal, n_variants))

    if effect_sizes is None:
        effect_sizes = [random.gauss(0, 0.5) for _ in range(len(causal_indices))]
    else:
        effect_sizes = list(effect_sizes[: len(causal_indices)])
        while len(effect_sizes) < len(causal_indices):
            effect_sizes.append(random.gauss(0, 0.5))

    # Generate phenotypes
    phenotypes: list[float] = []
    genetic_values: list[float] = []
    for i in range(n_samples):
        gv = sum(genotypes[i][ci] * es for ci, es in zip(causal_indices, effect_sizes))
        genetic_values.append(gv)

    gv_var = _std(genetic_values) ** 2 if len(genetic_values) > 1 else 1.0
    noise_var = gv_var  # h2 ~ 0.5
    heritability = gv_var / (gv_var + noise_var) if (gv_var + noise_var) > 0 else 0.0

    for gv in genetic_values:
        noise = random.gauss(0, math.sqrt(max(noise_var, 0.01)))
        phenotypes.append(gv + noise)

    logger.info(
        "Generated synthetic GWAS: %d samples, %d variants, %d causal, h2=%.3f",
        n_samples,
        n_variants,
        len(causal_indices),
        heritability,
    )

    return {
        "genotypes": genotypes,
        "phenotypes": phenotypes,
        "causal_variants": causal_indices,
        "true_effects": effect_sizes,
        "heritability": heritability,
        "maf": mafs,
    }


# ---------------------------------------------------------------------------
# Synthetic RNA-seq data
# ---------------------------------------------------------------------------


def generate_synthetic_expression(
    n_genes: int = 5000,
    n_samples: int = 50,
    n_de_genes: int = 200,
    fold_changes: list[float] | None = None,
    seed: int | None = None,
) -> dict:
    """Generate synthetic RNA-seq count data with known DE genes.

    Uses a negative binomial-like model to simulate count data with
    realistic mean-variance relationships.

    Args:
        n_genes: Total number of genes.
        n_samples: Total number of samples (split evenly between groups).
        n_de_genes: Number of differentially expressed genes.
        fold_changes: Optional list of log2 fold changes for DE genes.
        seed: Random seed.

    Returns:
        Dictionary with keys:
            - counts: 2D list ``(n_genes, n_samples)`` of count values.
            - groups: List of group labels (0 or 1) per sample.
            - de_genes: List of DE gene indices.
            - true_fold_changes: Log2 fold changes for DE genes.
    """
    _set_seed(seed)

    half = n_samples // 2
    groups = [0] * half + [1] * (n_samples - half)

    # Select DE genes
    de_indices = random.sample(range(n_genes), min(n_de_genes, n_genes))

    if fold_changes is None:
        fold_changes_list = [random.choice([-1, 1]) * random.uniform(0.5, 3.0) for _ in range(len(de_indices))]
    else:
        fold_changes_list = list(fold_changes[: len(de_indices)])
        while len(fold_changes_list) < len(de_indices):
            fold_changes_list.append(random.choice([-1, 1]) * random.uniform(0.5, 3.0))

    de_fc_map = dict(zip(de_indices, fold_changes_list))

    # Generate counts
    counts: list[list[int]] = []
    for g in range(n_genes):
        # Base mean expression (log-normal)
        base_mean = math.exp(random.gauss(3.0, 2.0))
        base_mean = max(1.0, base_mean)

        # Dispersion
        dispersion = 0.1 + random.random() * 0.5

        row: list[int] = []
        for s in range(n_samples):
            mean_val = base_mean
            if g in de_fc_map and groups[s] == 1:
                fc = de_fc_map[g]
                mean_val = base_mean * (2.0**fc)

            mean_val = max(0.1, mean_val)

            # Negative binomial as Gamma-Poisson
            # shape = 1/dispersion, scale = mean * dispersion
            shape = 1.0 / dispersion
            scale = mean_val * dispersion
            gamma_val = max(0.01, random.gammavariate(shape, scale))
            count = _poisson_sample(gamma_val)
            row.append(count)

        counts.append(row)

    logger.info(
        "Generated synthetic expression: %d genes, %d samples, %d DE genes",
        n_genes,
        n_samples,
        len(de_indices),
    )

    return {
        "counts": counts,
        "groups": groups,
        "de_genes": de_indices,
        "true_fold_changes": fold_changes_list,
    }


def _poisson_sample(lam: float) -> int:
    """Sample from Poisson distribution using Knuth's algorithm."""
    if lam <= 0:
        return 0
    if lam > 700:
        # For large lambda, use normal approximation
        return max(0, int(round(random.gauss(lam, math.sqrt(lam)))))

    l_val = math.exp(-lam)
    k = 0
    p = 1.0
    while p > l_val:
        k += 1
        p *= random.random()
    return k - 1


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------


def evaluate_benchmark(
    predictions: Any,
    truth: Any,
    task: str = "classification",
) -> dict:
    """Evaluate predictions against ground truth.

    Args:
        predictions: Predicted values/labels (list or array).
        truth: True values/labels (list or array).
        task: Task type: ``"classification"``, ``"regression"``, or
            ``"clustering"``.

    Returns:
        Dictionary with keys:
            - metrics: Dict of task-specific metric values.
            - confusion_matrix: For classification tasks (2D list).
            - summary: Human-readable summary.

    Raises:
        ValueError: If predictions and truth have different lengths.
    """
    pred_list = list(predictions)
    truth_list = list(truth)

    if len(pred_list) != len(truth_list):
        raise ValueError(f"predictions ({len(pred_list)}) and truth ({len(truth_list)}) " "must have same length")

    if task == "classification":
        return _evaluate_classification(pred_list, truth_list)
    elif task == "regression":
        return _evaluate_regression(pred_list, truth_list)
    elif task == "clustering":
        return _evaluate_clustering(pred_list, truth_list)
    else:
        raise ValueError(f"Unknown task '{task}'")


def _evaluate_classification(pred: list, truth: list) -> dict:
    """Evaluate classification predictions."""
    n = len(pred)
    correct = sum(1 for p, t in zip(pred, truth) if p == t)
    accuracy = correct / n if n > 0 else 0.0

    # Binary metrics
    classes = sorted(set(truth))
    if len(classes) == 2:
        pos = classes[1]
        tp = sum(1 for p, t in zip(pred, truth) if p == pos and t == pos)
        fp = sum(1 for p, t in zip(pred, truth) if p == pos and t != pos)
        fn = sum(1 for p, t in zip(pred, truth) if p != pos and t == pos)
        tn = sum(1 for p, t in zip(pred, truth) if p != pos and t != pos)

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

        confusion = [[tn, fp], [fn, tp]]
    else:
        precision = recall = f1 = 0.0
        confusion = []

    metrics = {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1_score": f1,
        "n_samples": n,
    }

    return {
        "metrics": metrics,
        "confusion_matrix": confusion,
        "summary": f"Classification: accuracy={accuracy:.3f}, F1={f1:.3f}",
    }


def _evaluate_regression(pred: list, truth: list) -> dict:
    """Evaluate regression predictions."""
    n = len(pred)
    residuals = [float(p) - float(t) for p, t in zip(pred, truth)]
    mse = sum(r**2 for r in residuals) / n if n > 0 else 0.0
    rmse = math.sqrt(mse)
    mae = sum(abs(r) for r in residuals) / n if n > 0 else 0.0

    truth_mean = _mean([float(t) for t in truth])
    ss_tot = sum((float(t) - truth_mean) ** 2 for t in truth)
    ss_res = sum(r**2 for r in residuals)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    metrics = {
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "r_squared": r_squared,
        "n_samples": n,
    }

    return {
        "metrics": metrics,
        "confusion_matrix": [],
        "summary": f"Regression: RMSE={rmse:.4f}, R2={r_squared:.3f}",
    }


def _evaluate_clustering(pred: list, truth: list) -> dict:
    """Evaluate clustering using Adjusted Rand Index approximation."""
    n = len(pred)

    # Contingency table
    pred_labels = sorted(set(pred))
    truth_labels = sorted(set(truth))
    contingency: dict[tuple, int] = {}
    for p, t in zip(pred, truth):
        key = (p, t)
        contingency[key] = contingency.get(key, 0) + 1

    # Rand Index
    tp_fp = 0
    for key, count in contingency.items():
        tp_fp += count * (count - 1) // 2

    # Row and column sums
    row_sums: dict[Any, int] = {}
    col_sums: dict[Any, int] = {}
    for (p, t), count in contingency.items():
        row_sums[p] = row_sums.get(p, 0) + count
        col_sums[t] = col_sums.get(t, 0) + count

    sum_row = sum(s * (s - 1) // 2 for s in row_sums.values())
    sum_col = sum(s * (s - 1) // 2 for s in col_sums.values())
    n_pairs = n * (n - 1) // 2

    expected = sum_row * sum_col / n_pairs if n_pairs > 0 else 0.0
    max_idx = (sum_row + sum_col) / 2.0

    if max_idx - expected == 0:
        ari = 0.0
    else:
        ari = (tp_fp - expected) / (max_idx - expected)

    metrics = {
        "adjusted_rand_index": ari,
        "n_predicted_clusters": len(pred_labels),
        "n_true_clusters": len(truth_labels),
        "n_samples": n,
    }

    return {
        "metrics": metrics,
        "confusion_matrix": [],
        "summary": f"Clustering: ARI={ari:.3f}, {len(pred_labels)} predicted / {len(truth_labels)} true clusters",
    }


# ---------------------------------------------------------------------------
# Benchmark suite
# ---------------------------------------------------------------------------


def benchmark_suite(
    methods: dict,
    dataset: dict,
    n_repeats: int = 5,
) -> dict:
    """Run multiple methods on a benchmark dataset and compare performance.

    Args:
        methods: Dict mapping method name to callable. Each callable takes
            ``(X, y)`` and returns predictions.
        dataset: Benchmark dataset from :func:`generate_benchmark_dataset`.
        n_repeats: Number of repeated evaluations (for variability estimate).

    Returns:
        Dictionary with keys:
            - results_per_method: Dict mapping method name to list of
              evaluation results.
            - rankings: List of (method_name, mean_score) tuples sorted by
              performance.
            - best_method: Name of the best-performing method.
            - summary: Comparison summary string.
    """
    task = dataset.get("metadata", {}).get("task", "classification")
    X = dataset["X"]
    y = dataset["y"]

    results_per_method: dict[str, list[dict]] = {}
    mean_scores: dict[str, float] = {}

    for method_name, method_fn in methods.items():
        method_results = []
        for rep in range(n_repeats):
            try:
                preds = method_fn(X, y)
                eval_result = evaluate_benchmark(preds, y, task=task)
                method_results.append(eval_result)
            except Exception as e:
                logger.warning(
                    "Method '%s' failed on repeat %d: %s",
                    method_name,
                    rep,
                    str(e),
                )
                method_results.append(
                    {
                        "metrics": {},
                        "confusion_matrix": [],
                        "summary": f"Failed: {e}",
                    }
                )

        results_per_method[method_name] = method_results

        # Extract primary metric
        primary_metric = _primary_metric_name(task)
        scores = [r["metrics"].get(primary_metric, 0.0) for r in method_results if r["metrics"]]
        mean_scores[method_name] = _mean(scores) if scores else 0.0

    # Rankings
    rankings = sorted(mean_scores.items(), key=lambda x: x[1], reverse=True)
    best = rankings[0][0] if rankings else "none"

    summary_lines = [f"Benchmark comparison ({task}, {n_repeats} repeats):"]
    for rank, (name, score) in enumerate(rankings, 1):
        summary_lines.append(f"  {rank}. {name}: {score:.4f}")

    logger.info("Benchmark suite complete: best=%s, %d methods compared", best, len(methods))

    return {
        "results_per_method": results_per_method,
        "rankings": rankings,
        "best_method": best,
        "summary": "\n".join(summary_lines),
    }


def _primary_metric_name(task: str) -> str:
    """Return primary metric name for a given task."""
    return {
        "classification": "accuracy",
        "regression": "r_squared",
        "clustering": "adjusted_rand_index",
        "de_genes": "accuracy",
    }.get(task, "accuracy")
