"""Pipeline stage definitions for long-read analysis workflows.

Contains the step-building functions for each pipeline type (QC, assembly,
methylation, SV, full). These are called by the LongReadOrchestrator methods
in orchestrator_core.py.
"""

from __future__ import annotations

import time
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any

from metainformant.core.utils.logging import get_logger

from .orchestrator_core import PipelineResult, PipelineStep, _extract_methylation_features

if TYPE_CHECKING:
    from .orchestrator_core import LongReadOrchestrator

logger = get_logger(__name__)


def _build_qc_steps(orch: LongReadOrchestrator, reads: Any) -> list[PipelineStep]:
    """Build QC pipeline steps.

    Steps:
        read_input -> filter_length -> filter_quality -> trim_adapters ->
        compute_metrics -> generate_qc_plots -> export_report
    """
    from ..quality.filtering import filter_by_length, filter_by_quality, trim_adapters
    from ..quality.metrics import calculate_n50, quality_score_distribution, read_length_stats

    params = orch.config.get("parameters", {})
    min_length = params.get("min_length", 1000)
    min_quality = params.get("min_quality", 7.0)
    do_trim = params.get("trim_adapters", True)
    max_search_length = params.get("max_search_length", 200)
    min_adapter_identity = params.get("min_adapter_identity", 0.75)

    steps: list[PipelineStep] = []

    # Step 1: Read input
    def _read_input(ctx: dict[str, Any]) -> list[dict[str, Any]]:
        return orch._load_reads_sync(reads)

    steps.append(PipelineStep(name="read_input", function=_read_input, params={}, depends_on=[]))

    # Step 2: Filter by length
    def _filter_length(ctx: dict[str, Any]) -> list[Any]:
        return filter_by_length(ctx["read_input"], min_length=min_length)

    steps.append(PipelineStep(name="filter_length", function=_filter_length, params={}, depends_on=["read_input"]))

    # Step 3: Filter by quality
    def _filter_quality(ctx: dict[str, Any]) -> list[Any]:
        return filter_by_quality(ctx["filter_length"], min_q=min_quality)

    steps.append(
        PipelineStep(name="filter_quality", function=_filter_quality, params={}, depends_on=["filter_length"])
    )

    # Step 4: Trim adapters (optional)
    if do_trim:

        def _trim_adapters(ctx: dict[str, Any]) -> list[Any]:
            return trim_adapters(
                ctx["filter_quality"],
                max_search_length=max_search_length,
                min_identity=min_adapter_identity,
            )

        steps.append(
            PipelineStep(name="trim_adapters", function=_trim_adapters, params={}, depends_on=["filter_quality"])
        )
        metric_dep = "trim_adapters"
    else:
        metric_dep = "filter_quality"

    # Step 5: Compute metrics
    def _compute_metrics(ctx: dict[str, Any]) -> dict[str, Any]:
        processed_reads = ctx[metric_dep]
        length_stats = read_length_stats(processed_reads)
        qual_dist = quality_score_distribution(processed_reads)
        return {
            "length_stats": length_stats,
            "quality_distribution": qual_dist,
            "total_reads_input": len(ctx["read_input"]),
            "reads_after_length_filter": len(ctx["filter_length"]),
            "reads_after_quality_filter": len(ctx["filter_quality"]),
            "reads_after_processing": len(processed_reads),
        }

    steps.append(PipelineStep(name="compute_metrics", function=_compute_metrics, params={}, depends_on=[metric_dep]))

    # Step 6: Generate QC plots
    def _generate_qc_plots(ctx: dict[str, Any]) -> dict[str, Path]:
        plots_dir = orch.output_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        generated: dict[str, Path] = {}

        processed_reads = ctx[metric_dep]
        try:
            from ..visualization.plots import plot_quality_vs_length, plot_read_length_histogram

            length_plot = plot_read_length_histogram(
                processed_reads,
                plots_dir / "read_length_histogram.png",
            )
            generated["read_length_histogram"] = length_plot

            quality_plot = plot_quality_vs_length(
                processed_reads,
                plots_dir / "quality_vs_length.png",
            )
            generated["quality_vs_length"] = quality_plot
        except ImportError:
            logger.warning("matplotlib not available; skipping QC plot generation")
        except Exception as exc:
            logger.warning("Plot generation failed: %s", exc)

        return generated

    steps.append(
        PipelineStep(name="generate_qc_plots", function=_generate_qc_plots, params={}, depends_on=["compute_metrics"])
    )

    # Step 7: Export report
    def _export_report(ctx: dict[str, Any]) -> Path:
        from .reporting import export_report, generate_qc_report

        metrics = ctx["compute_metrics"]
        interim_result = PipelineResult(
            pipeline_name="qc",
            steps=steps,
            success=True,
            total_duration=0.0,
            output_dir=orch.output_dir,
            summary=metrics,
        )
        qc_report = generate_qc_report(interim_result)
        return export_report(qc_report, orch.output_dir / "qc_report.json", format="json")

    steps.append(
        PipelineStep(name="export_report", function=_export_report, params={}, depends_on=["compute_metrics"])
    )

    return steps


def _build_assembly_steps(orch: LongReadOrchestrator, reads: Any) -> list[PipelineStep]:
    """Build assembly pipeline steps.

    Steps:
        filter_reads -> find_overlaps -> build_graph ->
        generate_consensus -> polish -> calculate_stats -> assembly_plots
    """
    from ..assembly.consensus import calculate_consensus_quality, generate_consensus, polish_consensus
    from ..assembly.overlap import compute_overlap_graph, filter_contained_reads, find_overlaps
    from ..quality.filtering import filter_by_length, filter_by_quality
    from ..quality.metrics import read_length_stats

    params = orch.config.get("parameters", {})
    min_overlap = params.get("min_overlap", 2000)
    k = params.get("k", 15)
    w = params.get("w", 10)
    polish_iterations = params.get("polish_iterations", 2)
    min_minimizer_matches = params.get("min_minimizer_matches", 3)
    max_overhang = params.get("max_overhang", 1000)
    min_read_length = params.get("min_read_length", 1000)
    min_read_quality = params.get("min_read_quality", 7.0)

    steps: list[PipelineStep] = []

    def _filter_reads(ctx: dict[str, Any]) -> list[dict[str, Any]]:
        raw_reads = orch._load_reads_sync(reads)
        length_filtered = filter_by_length(raw_reads, min_length=min_read_length)
        quality_filtered = filter_by_quality(length_filtered, min_q=min_read_quality)
        return quality_filtered  # type: ignore[return-value]

    steps.append(PipelineStep(name="filter_reads", function=_filter_reads, params={}, depends_on=[]))

    def _find_overlaps(ctx: dict[str, Any]) -> list[Any]:
        return find_overlaps(
            ctx["filter_reads"],
            min_overlap=min_overlap,
            k=k,
            w=w,
            min_minimizer_matches=min_minimizer_matches,
            max_overhang=max_overhang,
        )

    steps.append(
        PipelineStep(name="find_overlaps", function=_find_overlaps, params={}, depends_on=["filter_reads"])
    )

    def _build_graph(ctx: dict[str, Any]) -> dict[str, Any]:
        overlaps = ctx["find_overlaps"]
        filtered_overlaps = filter_contained_reads(overlaps)
        graph = compute_overlap_graph(filtered_overlaps)
        return {"graph": graph, "overlaps": filtered_overlaps}

    steps.append(PipelineStep(name="build_graph", function=_build_graph, params={}, depends_on=["find_overlaps"]))

    def _generate_consensus(ctx: dict[str, Any]) -> Any:
        filtered_reads = ctx["filter_reads"]
        if not filtered_reads:
            return None
        return generate_consensus(filtered_reads)

    steps.append(
        PipelineStep(name="generate_consensus", function=_generate_consensus, params={}, depends_on=["build_graph"])
    )

    def _polish(ctx: dict[str, Any]) -> Any:
        consensus_result = ctx["generate_consensus"]
        if consensus_result is None:
            return None
        return polish_consensus(
            consensus_result,
            ctx["filter_reads"],
            iterations=polish_iterations,
        )

    steps.append(PipelineStep(name="polish", function=_polish, params={}, depends_on=["generate_consensus"]))

    def _calculate_stats(ctx: dict[str, Any]) -> dict[str, Any]:
        polished = ctx["polish"]
        filtered_reads = ctx["filter_reads"]
        graph_data = ctx["build_graph"]

        stats: dict[str, Any] = {
            "num_input_reads": len(filtered_reads),
            "num_overlaps": len(graph_data["overlaps"]),
            "graph_nodes": graph_data["graph"].num_nodes,
            "graph_edges": graph_data["graph"].num_edges,
        }

        if polished is not None:
            stats["consensus_length"] = polished.length
            stats["consensus_num_reads"] = polished.num_reads
            stats["consensus_mean_quality"] = polished.mean_quality
            stats["consensus_mean_coverage"] = polished.mean_coverage

            qual = calculate_consensus_quality(polished, filtered_reads)
            stats["overall_quality"] = qual.get("overall_quality", 0.0)
            stats["overall_agreement"] = qual.get("overall_agreement", 0.0)

        input_stats = read_length_stats(filtered_reads)
        stats["input_n50"] = input_stats.n50
        stats["input_mean_length"] = input_stats.mean_length
        stats["input_total_bases"] = input_stats.total_bases

        return stats

    steps.append(PipelineStep(name="calculate_stats", function=_calculate_stats, params={}, depends_on=["polish"]))

    def _assembly_plots(ctx: dict[str, Any]) -> dict[str, Path]:
        plots_dir = orch.output_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        generated: dict[str, Path] = {}

        try:
            from ..visualization.plots import plot_read_length_histogram

            filtered_reads = ctx["filter_reads"]
            length_plot = plot_read_length_histogram(
                filtered_reads,
                plots_dir / "assembly_input_lengths.png",
                title="Assembly Input Read Lengths",
            )
            generated["input_length_histogram"] = length_plot
        except ImportError:
            logger.warning("matplotlib not available; skipping assembly plots")
        except Exception as exc:
            logger.warning("Assembly plot generation failed: %s", exc)

        return generated

    steps.append(
        PipelineStep(name="assembly_plots", function=_assembly_plots, params={}, depends_on=["calculate_stats"])
    )

    return steps


def _build_methylation_steps(orch: LongReadOrchestrator, reads: Any) -> list[PipelineStep]:
    """Build methylation pipeline steps.

    Steps:
        load_reads -> call_5mc -> call_6ma -> aggregate_regions ->
        differential_analysis -> methylation_plots -> methylation_report
    """
    from ..analysis.modified_bases import aggregate_methylation, call_5mc, call_6ma, differential_methylation

    params = orch.config.get("parameters", {})
    modification_types = params.get("modification_types", ["5mC", "6mA"])
    min_coverage = params.get("min_coverage", 5)
    significance_threshold = params.get("significance_threshold", 0.05)
    min_difference = params.get("min_difference", 0.2)
    methylation_threshold = params.get("methylation_threshold", 0.5)
    regions = params.get("regions", [])

    steps: list[PipelineStep] = []

    def _load_meth_reads(ctx: dict[str, Any]) -> list[dict[str, Any]]:
        return orch._load_reads_sync(reads)

    steps.append(PipelineStep(name="load_reads", function=_load_meth_reads, params={}, depends_on=[]))

    if "5mC" in modification_types:

        def _call_5mc(ctx: dict[str, Any]) -> list[Any]:
            loaded = ctx["load_reads"]
            features = _extract_methylation_features(loaded, mod_type="5mC")
            if not features:
                logger.info("No 5mC signal features found; skipping 5mC calling")
                return []
            return call_5mc(features)

        steps.append(PipelineStep(name="call_5mc", function=_call_5mc, params={}, depends_on=["load_reads"]))

    if "6mA" in modification_types:

        def _call_6ma(ctx: dict[str, Any]) -> list[Any]:
            loaded = ctx["load_reads"]
            features = _extract_methylation_features(loaded, mod_type="6mA")
            if not features:
                logger.info("No 6mA signal features found; skipping 6mA calling")
                return []
            return call_6ma(features)

        steps.append(PipelineStep(name="call_6ma", function=_call_6ma, params={}, depends_on=["load_reads"]))

    def _aggregate_regions(ctx: dict[str, Any]) -> list[Any]:
        all_calls = []
        if "call_5mc" in ctx and ctx["call_5mc"]:
            all_calls.extend(ctx["call_5mc"])
        if "call_6ma" in ctx and ctx["call_6ma"]:
            all_calls.extend(ctx["call_6ma"])

        if not all_calls or not regions:
            return []

        return aggregate_methylation(all_calls, regions, min_coverage=min_coverage)

    agg_deps = []
    if "5mC" in modification_types:
        agg_deps.append("call_5mc")
    if "6mA" in modification_types:
        agg_deps.append("call_6ma")

    steps.append(PipelineStep(name="aggregate_regions", function=_aggregate_regions, params={}, depends_on=agg_deps))

    def _differential_analysis(ctx: dict[str, Any]) -> dict[str, Any]:
        all_calls = []
        if "call_5mc" in ctx and ctx["call_5mc"]:
            all_calls.extend(ctx["call_5mc"])
        if "call_6ma" in ctx and ctx["call_6ma"]:
            all_calls.extend(ctx["call_6ma"])

        if not all_calls:
            return {"results": [], "num_tested": 0, "num_significant": 0}

        midpoint = len(all_calls) // 2
        if midpoint < 2:
            return {"results": [], "num_tested": 0, "num_significant": 0}

        sample1 = all_calls[:midpoint]
        sample2 = all_calls[midpoint:]

        diff_results = differential_methylation(
            sample1,
            sample2,
            regions=regions if regions else None,
            min_coverage=min_coverage,
            min_difference=min_difference,
            alpha=significance_threshold,
        )

        return {
            "results": diff_results,
            "num_tested": len(diff_results),
            "num_significant": sum(1 for r in diff_results if r.is_significant),
        }

    steps.append(
        PipelineStep(
            name="differential_analysis", function=_differential_analysis, params={}, depends_on=["aggregate_regions"]
        )
    )

    def _methylation_plots(ctx: dict[str, Any]) -> dict[str, Path]:
        plots_dir = orch.output_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        generated: dict[str, Path] = {}

        all_calls = []
        if "call_5mc" in ctx and ctx["call_5mc"]:
            all_calls.extend(ctx["call_5mc"])
        if "call_6ma" in ctx and ctx["call_6ma"]:
            all_calls.extend(ctx["call_6ma"])

        if not all_calls:
            return generated

        try:
            from ..visualization.plots import plot_methylation_track

            call_dicts = [
                {"position": c.position, "probability": c.probability, "chromosome": c.chromosome} for c in all_calls
            ]

            if regions:
                for i, region in enumerate(regions):
                    plot_path = plots_dir / f"methylation_region_{i}.png"
                    plot_methylation_track(call_dicts, region, plot_path)
                    generated[f"methylation_region_{i}"] = plot_path
            elif call_dicts:
                chroms = {c["chromosome"] for c in call_dicts if c["chromosome"]}
                positions = [c["position"] for c in call_dicts if c["position"] > 0]
                if chroms and positions:
                    auto_region = {
                        "chromosome": sorted(chroms)[0],
                        "start": min(positions),
                        "end": max(positions) + 1,
                    }
                    plot_path = plots_dir / "methylation_overview.png"
                    plot_methylation_track(call_dicts, auto_region, plot_path)
                    generated["methylation_overview"] = plot_path

        except ImportError:
            logger.warning("matplotlib not available; skipping methylation plots")
        except Exception as exc:
            logger.warning("Methylation plot generation failed: %s", exc)

        return generated

    steps.append(
        PipelineStep(
            name="methylation_plots", function=_methylation_plots, params={}, depends_on=["aggregate_regions"]
        )
    )

    def _methylation_report(ctx: dict[str, Any]) -> Path:
        from .reporting import export_report, generate_methylation_report

        interim_result = PipelineResult(
            pipeline_name="methylation",
            steps=steps,
            success=True,
            total_duration=0.0,
            output_dir=orch.output_dir,
            summary={
                "5mc_calls": len(ctx.get("call_5mc", [])),
                "6ma_calls": len(ctx.get("call_6ma", [])),
                "aggregated_regions": len(ctx.get("aggregate_regions", [])),
                "differential": ctx.get("differential_analysis", {}),
            },
        )
        report = generate_methylation_report(interim_result)
        return export_report(report, orch.output_dir / "methylation_report.json", format="json")

    report_deps = ["differential_analysis", "methylation_plots"]
    steps.append(
        PipelineStep(name="methylation_report", function=_methylation_report, params={}, depends_on=report_deps)
    )

    return steps


def _build_sv_steps(orch: LongReadOrchestrator, alignments: Any) -> list[PipelineStep]:
    """Build SV pipeline steps.

    Steps:
        load_alignments -> detect_svs -> detect_insertions ->
        detect_inversions -> phase_variants -> sv_summary -> sv_plots
    """
    from ..analysis.structural import (
        detect_insertions,
        detect_inversions,
        detect_sv_from_long_reads,
        phase_structural_variants,
    )

    params = orch.config.get("parameters", {})
    min_sv_size = params.get("min_sv_size", 50)
    min_support = params.get("min_support", 3)
    sv_types = params.get("sv_types", ["DEL", "INS", "INV", "DUP", "BND"])
    min_mapping_quality = params.get("min_mapping_quality", 20)

    steps: list[PipelineStep] = []

    def _load_sv_alignments(ctx: dict[str, Any]) -> list[dict[str, Any]]:
        return orch._load_alignments_sync(alignments)

    steps.append(PipelineStep(name="load_alignments", function=_load_sv_alignments, params={}, depends_on=[]))

    def _detect_svs(ctx: dict[str, Any]) -> list[Any]:
        return detect_sv_from_long_reads(
            ctx["load_alignments"],
            min_size=min_sv_size,
            min_support=min_support,
            min_mapping_quality=min_mapping_quality,
        )

    steps.append(PipelineStep(name="detect_svs", function=_detect_svs, params={}, depends_on=["load_alignments"]))

    if "INS" in sv_types:

        def _detect_insertions(ctx: dict[str, Any]) -> list[Any]:
            return detect_insertions(
                ctx["load_alignments"],
                min_size=min_sv_size,
                min_support=min_support,
            )

        steps.append(
            PipelineStep(
                name="detect_insertions", function=_detect_insertions, params={}, depends_on=["load_alignments"]
            )
        )

    if "INV" in sv_types:

        def _detect_inversions(ctx: dict[str, Any]) -> list[Any]:
            return detect_inversions(
                ctx["load_alignments"],
                min_size=min_sv_size,
                min_support=min_support,
            )

        steps.append(
            PipelineStep(
                name="detect_inversions", function=_detect_inversions, params={}, depends_on=["load_alignments"]
            )
        )

    def _phase_variants(ctx: dict[str, Any]) -> list[Any]:
        svs = ctx["detect_svs"]
        if not svs:
            return []

        haplotype_tags: dict[str, int] = {}
        for aln in ctx["load_alignments"]:
            read_name = ""
            hp_tag = 0
            if isinstance(aln, dict):
                read_name = aln.get("read_name", "")
                tags = aln.get("tags", {})
                hp_tag = tags.get("HP", 0) if isinstance(tags, dict) else 0
            elif hasattr(aln, "read_name"):
                read_name = aln.read_name
                if hasattr(aln, "tags") and isinstance(aln.tags, dict):
                    hp_tag = aln.tags.get("HP", 0)

            if read_name and hp_tag in (1, 2):
                haplotype_tags[read_name] = hp_tag

        return phase_structural_variants(svs, haplotype_tags)

    steps.append(PipelineStep(name="phase_variants", function=_phase_variants, params={}, depends_on=["detect_svs"]))

    def _sv_summary(ctx: dict[str, Any]) -> dict[str, Any]:
        all_svs = ctx["detect_svs"]
        phased_svs = ctx.get("phase_variants", [])
        insertions = ctx.get("detect_insertions", [])
        inversions = ctx.get("detect_inversions", [])

        type_counts: dict[str, int] = defaultdict(int)
        size_distributions: dict[str, list[int]] = defaultdict(list)
        for sv in all_svs:
            sv_type = sv.sv_type if hasattr(sv, "sv_type") else sv.get("sv_type", "UNKNOWN")
            sv_size = sv.size if hasattr(sv, "size") else sv.get("size", 0)
            type_counts[sv_type] += 1
            if sv_size > 0:
                size_distributions[sv_type].append(sv_size)

        size_stats: dict[str, dict[str, float]] = {}
        for sv_type, sizes in size_distributions.items():
            if sizes:
                sorted_sizes = sorted(sizes)
                n = len(sorted_sizes)
                size_stats[sv_type] = {
                    "count": n,
                    "min": sorted_sizes[0],
                    "max": sorted_sizes[-1],
                    "mean": sum(sorted_sizes) / n,
                    "median": sorted_sizes[n // 2],
                }

        phased_count = sum(
            1 for sv in phased_svs if (sv.haplotype if hasattr(sv, "haplotype") else sv.get("haplotype", 0)) > 0
        )

        return {
            "total_svs": len(all_svs),
            "type_counts": dict(type_counts),
            "size_stats": size_stats,
            "phased_svs": len(phased_svs),
            "phased_to_haplotype": phased_count,
            "specialized_insertions": len(insertions),
            "specialized_inversions": len(inversions),
        }

    summary_deps = ["detect_svs"]
    if "INS" in sv_types:
        summary_deps.append("detect_insertions")
    if "INV" in sv_types:
        summary_deps.append("detect_inversions")

    steps.append(PipelineStep(name="sv_summary", function=_sv_summary, params={}, depends_on=summary_deps))

    def _sv_plots(ctx: dict[str, Any]) -> dict[str, Path]:
        plots_dir = orch.output_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        generated: dict[str, Path] = {}

        all_svs = ctx["detect_svs"]
        if not all_svs:
            return generated

        try:
            from ..visualization.plots import plot_read_length_histogram

            sv_sizes = []
            for sv in all_svs:
                size = sv.size if hasattr(sv, "size") else sv.get("size", 0)
                if size > 0:
                    sv_sizes.append(size)

            if sv_sizes:
                size_plot = plot_read_length_histogram(
                    sv_sizes,
                    plots_dir / "sv_size_distribution.png",
                    title="Structural Variant Size Distribution",
                    log_scale=True,
                )
                generated["sv_size_distribution"] = size_plot

        except ImportError:
            logger.warning("matplotlib not available; skipping SV plots")
        except Exception as exc:
            logger.warning("SV plot generation failed: %s", exc)

        return generated

    steps.append(PipelineStep(name="sv_plots", function=_sv_plots, params={}, depends_on=["sv_summary"]))

    return steps


def _run_full_pipeline(orch: LongReadOrchestrator, reads: Any) -> PipelineResult:
    """Run the full pipeline: QC, Assembly, Methylation, SV calling.

    Executes all sub-pipelines in sequence, passing results forward.
    """
    from .pipelines import (
        get_assembly_pipeline_config,
        get_methylation_pipeline_config,
        get_qc_pipeline_config,
        get_sv_pipeline_config,
    )

    full_start = time.time()
    all_steps: list[PipelineStep] = []
    sub_results: dict[str, PipelineResult] = {}
    overall_success = True

    params = orch.config.get("parameters", {})

    # 1. QC Pipeline
    logger.info("Full pipeline: running QC sub-pipeline...")
    qc_config = get_qc_pipeline_config(
        min_length=params.get("min_length", 1000),
        min_quality=params.get("min_quality", 7.0),
        trim_adapters=params.get("trim_adapters", True),
    )

    # Import here to avoid circular dependency
    from .orchestrator_core import LongReadOrchestrator as _Orch

    qc_orchestrator = _Orch(qc_config, orch.output_dir / "qc")
    qc_result = qc_orchestrator.run_qc_pipeline(reads)
    sub_results["qc"] = qc_result
    all_steps.extend(qc_result.steps)
    if not qc_result.success:
        overall_success = False

    # Get filtered reads from QC for downstream steps
    filtered_reads = reads  # fallback
    for step in qc_result.steps:
        if step.name in ("trim_adapters", "filter_quality") and step.result is not None:
            filtered_reads = step.result

    # 2. Assembly Pipeline
    logger.info("Full pipeline: running assembly sub-pipeline...")
    assembly_config = get_assembly_pipeline_config(
        min_overlap=params.get("min_overlap", 2000),
        k=params.get("k", 15),
        w=params.get("w", 10),
        polish_iterations=params.get("polish_iterations", 2),
    )
    assembly_orchestrator = _Orch(assembly_config, orch.output_dir / "assembly")
    assembly_result = assembly_orchestrator.run_assembly_pipeline(filtered_reads)
    sub_results["assembly"] = assembly_result
    all_steps.extend(assembly_result.steps)
    if not assembly_result.success:
        overall_success = False

    # 3. Methylation Pipeline
    logger.info("Full pipeline: running methylation sub-pipeline...")
    methylation_config = get_methylation_pipeline_config(
        modification_types=params.get("modification_types", ["5mC", "6mA"]),
        min_coverage=params.get("min_coverage", 5),
        significance_threshold=params.get("significance_threshold", 0.05),
    )
    methylation_orchestrator = _Orch(methylation_config, orch.output_dir / "methylation")
    methylation_result = methylation_orchestrator.run_methylation_pipeline(filtered_reads)
    sub_results["methylation"] = methylation_result
    all_steps.extend(methylation_result.steps)
    if not methylation_result.success:
        overall_success = False

    # 4. SV Pipeline
    logger.info("Full pipeline: running SV sub-pipeline...")
    sv_config = get_sv_pipeline_config(
        min_sv_size=params.get("min_sv_size", 50),
        min_support=params.get("min_support", 3),
    )
    sv_orchestrator = _Orch(sv_config, orch.output_dir / "sv")
    sv_result = sv_orchestrator.run_sv_pipeline(reads)  # use original reads for SV
    sub_results["sv"] = sv_result
    all_steps.extend(sv_result.steps)
    if not sv_result.success:
        overall_success = False

    total_duration = time.time() - full_start

    summary = {
        "sub_pipelines": {
            name: {
                "success": result.success,
                "duration": result.total_duration,
                "steps_completed": sum(1 for s in result.steps if s.status == "completed"),
                "steps_failed": sum(1 for s in result.steps if s.status == "failed"),
            }
            for name, result in sub_results.items()
        },
        "total_duration": total_duration,
    }

    logger.info(
        "Full pipeline completed in %.1fs: success=%s",
        total_duration,
        overall_success,
    )

    return PipelineResult(
        pipeline_name="full",
        steps=all_steps,
        success=overall_success,
        total_duration=total_duration,
        output_dir=orch.output_dir,
        summary=summary,
    )
