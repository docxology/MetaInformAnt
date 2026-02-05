"""Pipeline orchestration engine for long-read analysis workflows.

Provides the LongReadOrchestrator class that executes end-to-end analysis
pipelines by chaining together functions from the longread subpackages
(io, quality, analysis, assembly, visualization). Handles dependency
resolution, step execution with timing, error handling, and result
aggregation.

Each pipeline is defined as a directed acyclic graph (DAG) of PipelineStep
objects. The orchestrator performs topological sorting to determine execution
order and propagates results through a shared context dictionary.
"""

from __future__ import annotations

import json
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


@dataclass
class PipelineStep:
    """A single step in an analysis pipeline.

    Attributes:
        name: Unique step identifier.
        function: Callable to execute for this step.
        params: Parameters to pass to the function.
        depends_on: Names of steps that must complete before this one.
        status: Current execution status.
        result: Return value from the function after execution.
        duration_seconds: Wall-clock execution time.
        error: Error message if the step failed.
    """

    name: str
    function: Callable[..., Any]
    params: dict[str, Any] = field(default_factory=dict)
    depends_on: list[str] = field(default_factory=list)
    status: str = "pending"  # pending | running | completed | failed | skipped
    result: Any = None
    duration_seconds: float = 0.0
    error: str = ""


@dataclass
class PipelineResult:
    """Aggregated result from a complete pipeline execution.

    Attributes:
        pipeline_name: Name of the executed pipeline.
        steps: All pipeline steps with their results and statuses.
        success: Whether all non-skipped steps completed successfully.
        total_duration: Total wall-clock duration of the pipeline.
        output_dir: Directory where outputs were written.
        summary: Pipeline-specific summary metrics.
    """

    pipeline_name: str
    steps: list[PipelineStep]
    success: bool
    total_duration: float
    output_dir: Path
    summary: dict[str, Any] = field(default_factory=dict)


class LongReadOrchestrator:
    """Orchestrates end-to-end long-read analysis pipelines.

    Ties together the io, quality, analysis, assembly, and visualization
    subpackages into cohesive pipelines. Each pipeline is defined as a
    sequence of dependent steps that are executed in topological order.

    Usage:
        config = get_qc_pipeline_config(min_length=1000, min_quality=7.0)
        orchestrator = LongReadOrchestrator(config, output_dir="output/longread/qc")
        result = orchestrator.run_qc_pipeline(reads)
    """

    def __init__(self, config: dict[str, Any], output_dir: Path | str) -> None:
        """Initialize the orchestrator.

        Args:
            config: Pipeline configuration dictionary containing parameters
                and step definitions. Typically from get_*_pipeline_config().
            output_dir: Directory for pipeline output files.
        """
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._context: dict[str, Any] = {}
        logger.info(
            "Initialized LongReadOrchestrator: pipeline=%s, output_dir=%s",
            config.get("pipeline_name", "unknown"),
            self.output_dir,
        )

    def run_qc_pipeline(self, reads: list[dict[str, Any]] | Path) -> PipelineResult:
        """Run the QC pipeline: load, filter, metrics, QC report, plots.

        Steps:
            read_input -> filter_length -> filter_quality -> trim_adapters ->
            compute_metrics -> generate_qc_plots -> export_report

        Args:
            reads: List of read dictionaries (with 'sequence' and optionally
                'quality_string', 'read_id' keys), or a Path to a FAST5/BAM file.

        Returns:
            PipelineResult with QC metrics and report paths.
        """
        from ..quality.filtering import filter_by_length, filter_by_quality, trim_adapters
        from ..quality.metrics import calculate_n50, quality_score_distribution, read_length_stats

        params = self.config.get("parameters", {})
        min_length = params.get("min_length", 1000)
        min_quality = params.get("min_quality", 7.0)
        do_trim = params.get("trim_adapters", True)
        max_search_length = params.get("max_search_length", 200)
        min_adapter_identity = params.get("min_adapter_identity", 0.75)

        steps: list[PipelineStep] = []

        # Step 1: Read input
        def _read_input(ctx: dict[str, Any]) -> list[dict[str, Any]]:
            return self._load_reads_sync(reads)

        steps.append(
            PipelineStep(
                name="read_input",
                function=_read_input,
                params={},
                depends_on=[],
            )
        )

        # Step 2: Filter by length
        def _filter_length(ctx: dict[str, Any]) -> list[Any]:
            return filter_by_length(ctx["read_input"], min_length=min_length)

        steps.append(
            PipelineStep(
                name="filter_length",
                function=_filter_length,
                params={},
                depends_on=["read_input"],
            )
        )

        # Step 3: Filter by quality
        def _filter_quality(ctx: dict[str, Any]) -> list[Any]:
            return filter_by_quality(ctx["filter_length"], min_q=min_quality)

        steps.append(
            PipelineStep(
                name="filter_quality",
                function=_filter_quality,
                params={},
                depends_on=["filter_length"],
            )
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
                PipelineStep(
                    name="trim_adapters",
                    function=_trim_adapters,
                    params={},
                    depends_on=["filter_quality"],
                )
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

        steps.append(
            PipelineStep(
                name="compute_metrics",
                function=_compute_metrics,
                params={},
                depends_on=[metric_dep],
            )
        )

        # Step 6: Generate QC plots
        def _generate_qc_plots(ctx: dict[str, Any]) -> dict[str, Path]:
            plots_dir = self.output_dir / "plots"
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
            PipelineStep(
                name="generate_qc_plots",
                function=_generate_qc_plots,
                params={},
                depends_on=["compute_metrics"],
            )
        )

        # Step 7: Export report
        def _export_report(ctx: dict[str, Any]) -> Path:
            from .reporting import export_report, generate_qc_report

            metrics = ctx["compute_metrics"]
            # Build a lightweight PipelineResult for the report generator
            interim_result = PipelineResult(
                pipeline_name="qc",
                steps=steps,
                success=True,
                total_duration=0.0,
                output_dir=self.output_dir,
                summary=metrics,
            )
            qc_report = generate_qc_report(interim_result)
            return export_report(qc_report, self.output_dir / "qc_report.json", format="json")

        steps.append(
            PipelineStep(
                name="export_report",
                function=_export_report,
                params={},
                depends_on=["compute_metrics"],
            )
        )

        return self._run_steps("qc", steps)

    def run_assembly_pipeline(self, reads: list[dict[str, Any]] | Path) -> PipelineResult:
        """Run the assembly pipeline: QC, overlap, consensus, polish, stats.

        Steps:
            filter_reads -> find_overlaps -> build_graph ->
            generate_consensus -> polish -> calculate_stats -> assembly_plots

        Args:
            reads: List of read dictionaries or Path to a FAST5/BAM file.

        Returns:
            PipelineResult with assembly contigs and statistics.
        """
        from ..assembly.consensus import calculate_consensus_quality, generate_consensus, polish_consensus
        from ..assembly.overlap import compute_overlap_graph, filter_contained_reads, find_overlaps
        from ..quality.filtering import filter_by_length, filter_by_quality
        from ..quality.metrics import read_length_stats

        params = self.config.get("parameters", {})
        min_overlap = params.get("min_overlap", 2000)
        k = params.get("k", 15)
        w = params.get("w", 10)
        polish_iterations = params.get("polish_iterations", 2)
        min_minimizer_matches = params.get("min_minimizer_matches", 3)
        max_overhang = params.get("max_overhang", 1000)
        min_read_length = params.get("min_read_length", 1000)
        min_read_quality = params.get("min_read_quality", 7.0)

        steps: list[PipelineStep] = []

        # Step 1: Filter reads
        def _filter_reads(ctx: dict[str, Any]) -> list[dict[str, Any]]:
            raw_reads = self._load_reads_sync(reads)
            length_filtered = filter_by_length(raw_reads, min_length=min_read_length)
            quality_filtered = filter_by_quality(length_filtered, min_q=min_read_quality)
            return quality_filtered  # type: ignore[return-value]

        steps.append(
            PipelineStep(
                name="filter_reads",
                function=_filter_reads,
                params={},
                depends_on=[],
            )
        )

        # Step 2: Find overlaps
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
            PipelineStep(
                name="find_overlaps",
                function=_find_overlaps,
                params={},
                depends_on=["filter_reads"],
            )
        )

        # Step 3: Build overlap graph
        def _build_graph(ctx: dict[str, Any]) -> dict[str, Any]:
            overlaps = ctx["find_overlaps"]
            filtered_overlaps = filter_contained_reads(overlaps)
            graph = compute_overlap_graph(filtered_overlaps)
            return {"graph": graph, "overlaps": filtered_overlaps}

        steps.append(
            PipelineStep(
                name="build_graph",
                function=_build_graph,
                params={},
                depends_on=["find_overlaps"],
            )
        )

        # Step 4: Generate consensus
        def _generate_consensus(ctx: dict[str, Any]) -> Any:
            filtered_reads = ctx["filter_reads"]
            if not filtered_reads:
                return None
            return generate_consensus(filtered_reads)

        steps.append(
            PipelineStep(
                name="generate_consensus",
                function=_generate_consensus,
                params={},
                depends_on=["build_graph"],
            )
        )

        # Step 5: Polish consensus
        def _polish(ctx: dict[str, Any]) -> Any:
            consensus_result = ctx["generate_consensus"]
            if consensus_result is None:
                return None
            return polish_consensus(
                consensus_result,
                ctx["filter_reads"],
                iterations=polish_iterations,
            )

        steps.append(
            PipelineStep(
                name="polish",
                function=_polish,
                params={},
                depends_on=["generate_consensus"],
            )
        )

        # Step 6: Calculate stats
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

                # Quality assessment
                qual = calculate_consensus_quality(polished, filtered_reads)
                stats["overall_quality"] = qual.get("overall_quality", 0.0)
                stats["overall_agreement"] = qual.get("overall_agreement", 0.0)

            # Input read stats
            input_stats = read_length_stats(filtered_reads)
            stats["input_n50"] = input_stats.n50
            stats["input_mean_length"] = input_stats.mean_length
            stats["input_total_bases"] = input_stats.total_bases

            return stats

        steps.append(
            PipelineStep(
                name="calculate_stats",
                function=_calculate_stats,
                params={},
                depends_on=["polish"],
            )
        )

        # Step 7: Assembly plots
        def _assembly_plots(ctx: dict[str, Any]) -> dict[str, Path]:
            plots_dir = self.output_dir / "plots"
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
            PipelineStep(
                name="assembly_plots",
                function=_assembly_plots,
                params={},
                depends_on=["calculate_stats"],
            )
        )

        return self._run_steps("assembly", steps)

    def run_methylation_pipeline(self, reads: list[dict[str, Any]] | Path) -> PipelineResult:
        """Run the methylation pipeline: load, call mods, aggregate, diff, plots.

        Steps:
            load_reads -> call_5mc -> call_6ma -> aggregate_regions ->
            differential_analysis -> methylation_plots -> methylation_report

        Args:
            reads: List of read dictionaries with signal features or
                methylation tags, or a Path to a BAM file.

        Returns:
            PipelineResult with methylation calls and analysis results.
        """
        from ..analysis.modified_bases import aggregate_methylation, call_5mc, call_6ma, differential_methylation

        params = self.config.get("parameters", {})
        modification_types = params.get("modification_types", ["5mC", "6mA"])
        min_coverage = params.get("min_coverage", 5)
        significance_threshold = params.get("significance_threshold", 0.05)
        min_difference = params.get("min_difference", 0.2)
        methylation_threshold = params.get("methylation_threshold", 0.5)
        regions = params.get("regions", [])

        steps: list[PipelineStep] = []

        # Step 1: Load reads
        def _load_meth_reads(ctx: dict[str, Any]) -> list[dict[str, Any]]:
            return self._load_reads_sync(reads)

        steps.append(
            PipelineStep(
                name="load_reads",
                function=_load_meth_reads,
                params={},
                depends_on=[],
            )
        )

        # Step 2: Call 5mC
        def _call_5mc(ctx: dict[str, Any]) -> list[Any]:
            loaded = ctx["load_reads"]
            # Extract signal features from reads for 5mC calling
            features = _extract_methylation_features(loaded, mod_type="5mC")
            if not features:
                logger.info("No 5mC signal features found; skipping 5mC calling")
                return []
            return call_5mc(features)

        if "5mC" in modification_types:
            steps.append(
                PipelineStep(
                    name="call_5mc",
                    function=_call_5mc,
                    params={},
                    depends_on=["load_reads"],
                )
            )

        # Step 3: Call 6mA
        def _call_6ma(ctx: dict[str, Any]) -> list[Any]:
            loaded = ctx["load_reads"]
            features = _extract_methylation_features(loaded, mod_type="6mA")
            if not features:
                logger.info("No 6mA signal features found; skipping 6mA calling")
                return []
            return call_6ma(features)

        if "6mA" in modification_types:
            steps.append(
                PipelineStep(
                    name="call_6ma",
                    function=_call_6ma,
                    params={},
                    depends_on=["load_reads"],
                )
            )

        # Step 4: Aggregate regions
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

        steps.append(
            PipelineStep(
                name="aggregate_regions",
                function=_aggregate_regions,
                params={},
                depends_on=agg_deps,
            )
        )

        # Step 5: Differential analysis
        def _differential_analysis(ctx: dict[str, Any]) -> dict[str, Any]:
            all_calls = []
            if "call_5mc" in ctx and ctx["call_5mc"]:
                all_calls.extend(ctx["call_5mc"])
            if "call_6ma" in ctx and ctx["call_6ma"]:
                all_calls.extend(ctx["call_6ma"])

            if not all_calls:
                return {"results": [], "num_tested": 0, "num_significant": 0}

            # Split calls into two halves for differential analysis demonstration
            # In real usage, two separate samples would be provided via config
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
                name="differential_analysis",
                function=_differential_analysis,
                params={},
                depends_on=["aggregate_regions"],
            )
        )

        # Step 6: Methylation plots
        def _methylation_plots(ctx: dict[str, Any]) -> dict[str, Path]:
            plots_dir = self.output_dir / "plots"
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

                # Plot methylation for each region if defined
                call_dicts = [
                    {"position": c.position, "probability": c.probability, "chromosome": c.chromosome}
                    for c in all_calls
                ]

                if regions:
                    for i, region in enumerate(regions):
                        plot_path = plots_dir / f"methylation_region_{i}.png"
                        plot_methylation_track(call_dicts, region, plot_path)
                        generated[f"methylation_region_{i}"] = plot_path
                elif call_dicts:
                    # Auto-detect region from calls
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
                name="methylation_plots",
                function=_methylation_plots,
                params={},
                depends_on=["aggregate_regions"],
            )
        )

        # Step 7: Methylation report
        def _methylation_report(ctx: dict[str, Any]) -> Path:
            from .reporting import export_report, generate_methylation_report

            interim_result = PipelineResult(
                pipeline_name="methylation",
                steps=steps,
                success=True,
                total_duration=0.0,
                output_dir=self.output_dir,
                summary={
                    "5mc_calls": len(ctx.get("call_5mc", [])),
                    "6ma_calls": len(ctx.get("call_6ma", [])),
                    "aggregated_regions": len(ctx.get("aggregate_regions", [])),
                    "differential": ctx.get("differential_analysis", {}),
                },
            )
            report = generate_methylation_report(interim_result)
            return export_report(report, self.output_dir / "methylation_report.json", format="json")

        report_deps = ["differential_analysis", "methylation_plots"]
        steps.append(
            PipelineStep(
                name="methylation_report",
                function=_methylation_report,
                params={},
                depends_on=report_deps,
            )
        )

        return self._run_steps("methylation", steps)

    def run_sv_pipeline(self, alignments: list[dict[str, Any]] | Path) -> PipelineResult:
        """Run the SV pipeline: load, detect SVs, phase, annotate, plots.

        Steps:
            load_alignments -> detect_svs -> detect_insertions ->
            detect_inversions -> phase_variants -> sv_summary -> sv_plots

        Args:
            alignments: List of alignment dictionaries or Path to a BAM file.

        Returns:
            PipelineResult with structural variant calls and statistics.
        """
        from ..analysis.structural import (
            detect_insertions,
            detect_inversions,
            detect_sv_from_long_reads,
            phase_structural_variants,
        )

        params = self.config.get("parameters", {})
        min_sv_size = params.get("min_sv_size", 50)
        min_support = params.get("min_support", 3)
        sv_types = params.get("sv_types", ["DEL", "INS", "INV", "DUP", "BND"])
        min_mapping_quality = params.get("min_mapping_quality", 20)

        steps: list[PipelineStep] = []

        # Step 1: Load alignments
        def _load_sv_alignments(ctx: dict[str, Any]) -> list[dict[str, Any]]:
            return self._load_alignments_sync(alignments)

        steps.append(
            PipelineStep(
                name="load_alignments",
                function=_load_sv_alignments,
                params={},
                depends_on=[],
            )
        )

        # Step 2: Detect all SVs
        def _detect_svs(ctx: dict[str, Any]) -> list[Any]:
            return detect_sv_from_long_reads(
                ctx["load_alignments"],
                min_size=min_sv_size,
                min_support=min_support,
                min_mapping_quality=min_mapping_quality,
            )

        steps.append(
            PipelineStep(
                name="detect_svs",
                function=_detect_svs,
                params={},
                depends_on=["load_alignments"],
            )
        )

        # Step 3: Detect insertions
        if "INS" in sv_types:

            def _detect_insertions(ctx: dict[str, Any]) -> list[Any]:
                return detect_insertions(
                    ctx["load_alignments"],
                    min_size=min_sv_size,
                    min_support=min_support,
                )

            steps.append(
                PipelineStep(
                    name="detect_insertions",
                    function=_detect_insertions,
                    params={},
                    depends_on=["load_alignments"],
                )
            )

        # Step 4: Detect inversions
        if "INV" in sv_types:

            def _detect_inversions(ctx: dict[str, Any]) -> list[Any]:
                return detect_inversions(
                    ctx["load_alignments"],
                    min_size=min_sv_size,
                    min_support=min_support,
                )

            steps.append(
                PipelineStep(
                    name="detect_inversions",
                    function=_detect_inversions,
                    params={},
                    depends_on=["load_alignments"],
                )
            )

        # Step 5: Phase variants
        def _phase_variants(ctx: dict[str, Any]) -> list[Any]:
            svs = ctx["detect_svs"]
            if not svs:
                return []

            # Build haplotype tags from alignment HP tags
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

        steps.append(
            PipelineStep(
                name="phase_variants",
                function=_phase_variants,
                params={},
                depends_on=["detect_svs"],
            )
        )

        # Step 6: SV summary
        def _sv_summary(ctx: dict[str, Any]) -> dict[str, Any]:
            all_svs = ctx["detect_svs"]
            phased_svs = ctx.get("phase_variants", [])
            insertions = ctx.get("detect_insertions", [])
            inversions = ctx.get("detect_inversions", [])

            # Count by type
            type_counts: dict[str, int] = defaultdict(int)
            size_distributions: dict[str, list[int]] = defaultdict(list)
            for sv in all_svs:
                sv_type = sv.sv_type if hasattr(sv, "sv_type") else sv.get("sv_type", "UNKNOWN")
                sv_size = sv.size if hasattr(sv, "size") else sv.get("size", 0)
                type_counts[sv_type] += 1
                if sv_size > 0:
                    size_distributions[sv_type].append(sv_size)

            # Compute size stats per type
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

        steps.append(
            PipelineStep(
                name="sv_summary",
                function=_sv_summary,
                params={},
                depends_on=summary_deps,
            )
        )

        # Step 7: SV plots
        def _sv_plots(ctx: dict[str, Any]) -> dict[str, Path]:
            plots_dir = self.output_dir / "plots"
            plots_dir.mkdir(parents=True, exist_ok=True)
            generated: dict[str, Path] = {}

            all_svs = ctx["detect_svs"]
            if not all_svs:
                return generated

            try:
                from ..visualization.plots import plot_read_length_histogram

                # Plot SV size distribution using the length histogram
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

        steps.append(
            PipelineStep(
                name="sv_plots",
                function=_sv_plots,
                params={},
                depends_on=["sv_summary"],
            )
        )

        return self._run_steps("sv", steps)

    def run_full_pipeline(self, reads: list[dict[str, Any]] | Path) -> PipelineResult:
        """Run the full pipeline: QC, Assembly, Methylation, SV calling.

        Executes all sub-pipelines in sequence, passing results forward.
        The full pipeline first runs QC to filter reads, then uses the
        filtered reads for assembly and methylation analysis. SV calling
        runs on the original alignment data.

        Args:
            reads: List of read dictionaries or Path to input file.

        Returns:
            PipelineResult with combined results from all sub-pipelines.
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

        params = self.config.get("parameters", {})

        # 1. QC Pipeline
        logger.info("Full pipeline: running QC sub-pipeline...")
        qc_config = get_qc_pipeline_config(
            min_length=params.get("min_length", 1000),
            min_quality=params.get("min_quality", 7.0),
            trim_adapters=params.get("trim_adapters", True),
        )
        qc_orchestrator = LongReadOrchestrator(qc_config, self.output_dir / "qc")
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
        assembly_orchestrator = LongReadOrchestrator(assembly_config, self.output_dir / "assembly")
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
        methylation_orchestrator = LongReadOrchestrator(methylation_config, self.output_dir / "methylation")
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
        sv_orchestrator = LongReadOrchestrator(sv_config, self.output_dir / "sv")
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
            output_dir=self.output_dir,
            summary=summary,
        )

    def run_pipeline(self, pipeline_name: str, input_data: Any) -> PipelineResult:
        """Run a named pipeline by dispatching to the appropriate method.

        Args:
            pipeline_name: One of 'qc', 'assembly', 'methylation', 'sv', 'full'.
            input_data: Input data (reads or alignments depending on pipeline).

        Returns:
            PipelineResult from the executed pipeline.

        Raises:
            ValueError: If pipeline_name is not recognized.
        """
        dispatch: dict[str, Callable[..., PipelineResult]] = {
            "qc": self.run_qc_pipeline,
            "assembly": self.run_assembly_pipeline,
            "methylation": self.run_methylation_pipeline,
            "sv": self.run_sv_pipeline,
            "full": self.run_full_pipeline,
        }

        if pipeline_name not in dispatch:
            raise ValueError(f"Unknown pipeline '{pipeline_name}'. Valid pipelines: {sorted(dispatch.keys())}")

        logger.info("Running pipeline: %s", pipeline_name)
        return dispatch[pipeline_name](input_data)

    # --- Internal methods ---

    def _run_steps(self, pipeline_name: str, steps: list[PipelineStep]) -> PipelineResult:
        """Execute a list of pipeline steps in dependency order.

        Args:
            pipeline_name: Name for the pipeline result.
            steps: List of PipelineStep objects to execute.

        Returns:
            PipelineResult with all step outcomes.
        """
        start_time = time.time()
        context: dict[str, Any] = {}

        # Topological sort
        ordered_steps = self._resolve_dependencies(steps)

        all_success = True

        for step in ordered_steps:
            # Check if dependencies completed
            deps_met = all(context.get(dep + "_status") == "completed" for dep in step.depends_on)

            if not deps_met:
                step.status = "skipped"
                step.error = "Dependencies not met"
                logger.warning("Skipping step '%s': dependencies not met", step.name)
                context[step.name + "_status"] = "skipped"
                continue

            # Execute the step
            result = self._execute_step(step, context)

            if step.status == "completed":
                context[step.name] = result
                context[step.name + "_status"] = "completed"
            else:
                context[step.name + "_status"] = "failed"
                all_success = False

        total_duration = time.time() - start_time

        # Build summary from completed steps
        summary: dict[str, Any] = {}
        for step in ordered_steps:
            if step.status == "completed" and isinstance(step.result, dict):
                summary[step.name] = step.result

        # Export pipeline summary to output dir
        summary_path = self.output_dir / f"{pipeline_name}_summary.json"
        try:
            _write_json_safe(summary, summary_path)
        except Exception as exc:
            logger.warning("Could not write pipeline summary: %s", exc)

        logger.info(
            "Pipeline '%s' completed in %.1fs: %d/%d steps succeeded",
            pipeline_name,
            total_duration,
            sum(1 for s in ordered_steps if s.status == "completed"),
            len(ordered_steps),
        )

        return PipelineResult(
            pipeline_name=pipeline_name,
            steps=ordered_steps,
            success=all_success,
            total_duration=total_duration,
            output_dir=self.output_dir,
            summary=summary,
        )

    def _execute_step(self, step: PipelineStep, context: dict[str, Any]) -> Any:
        """Execute a single pipeline step with timing and error handling.

        Args:
            step: The PipelineStep to execute.
            context: Shared context dictionary with results from prior steps.

        Returns:
            The return value of the step function.
        """
        logger.info("Executing step: %s", step.name)
        step.status = "running"
        start = time.time()

        try:
            result = step.function(context)
            step.status = "completed"
            step.result = result
            step.duration_seconds = time.time() - start
            logger.info(
                "Step '%s' completed in %.2fs",
                step.name,
                step.duration_seconds,
            )
            return result
        except Exception as exc:
            step.status = "failed"
            step.error = str(exc)
            step.duration_seconds = time.time() - start
            logger.error(
                "Step '%s' failed after %.2fs: %s",
                step.name,
                step.duration_seconds,
                exc,
            )
            return None

    def _resolve_dependencies(self, steps: list[PipelineStep]) -> list[PipelineStep]:
        """Topological sort of steps based on dependencies.

        Uses Kahn's algorithm for topological ordering. Steps with no
        dependencies are processed first.

        Args:
            steps: List of PipelineStep objects.

        Returns:
            List of PipelineStep objects in execution order.

        Raises:
            ValueError: If there is a circular dependency.
        """
        step_map: dict[str, PipelineStep] = {s.name: s for s in steps}
        in_degree: dict[str, int] = {s.name: 0 for s in steps}

        for step in steps:
            for dep in step.depends_on:
                if dep in step_map:
                    in_degree[step.name] += 1

        # Initialize queue with steps having no dependencies
        queue: list[str] = [name for name, degree in in_degree.items() if degree == 0]
        ordered: list[PipelineStep] = []

        while queue:
            # Sort queue for deterministic order
            queue.sort()
            current_name = queue.pop(0)
            ordered.append(step_map[current_name])

            # Reduce in-degree for dependent steps
            for step in steps:
                if current_name in step.depends_on:
                    in_degree[step.name] -= 1
                    if in_degree[step.name] == 0:
                        queue.append(step.name)

        if len(ordered) != len(steps):
            processed_names = {s.name for s in ordered}
            unprocessed = [s.name for s in steps if s.name not in processed_names]
            raise ValueError(f"Circular dependency detected in steps: {unprocessed}")

        return ordered

    def _load_reads(self, context: dict[str, Any]) -> list[dict[str, Any]]:
        """Load reads from various input sources.

        Handles list[dict], Path to FAST5, and Path to BAM.

        Args:
            context: Step context (not used, reads come from params via closure).

        Returns:
            List of read dictionaries.
        """
        reads = context.get("reads", context)
        return self._load_reads_sync(reads)

    def _load_reads_sync(self, reads: Any) -> list[dict[str, Any]]:
        """Synchronous read loading from various sources.

        Args:
            reads: List of dicts, Path to file, or other iterable.

        Returns:
            List of read dictionaries.
        """
        if isinstance(reads, (list, tuple)):
            # Convert to list of dicts if not already
            result: list[dict[str, Any]] = []
            for r in reads:
                if isinstance(r, dict):
                    result.append(r)
                elif hasattr(r, "sequence"):
                    # Convert dataclass-like objects to dicts
                    d: dict[str, Any] = {
                        "read_id": getattr(r, "read_id", getattr(r, "read_name", "")),
                        "sequence": r.sequence,
                    }
                    if hasattr(r, "quality_string") and r.quality_string:
                        d["quality_string"] = r.quality_string
                    result.append(d)
            logger.info("Loaded %d reads from input list", len(result))
            return result

        if isinstance(reads, Path) or (isinstance(reads, str) and Path(reads).exists()):
            reads_path = Path(reads)
            suffix = reads_path.suffix.lower()

            if suffix in (".fast5", ".pod5"):
                from ..io.fast5 import read_fast5

                fast5_reads = read_fast5(reads_path)
                return [
                    {
                        "read_id": r.read_id,
                        "sequence": r.sequence or "",
                        "quality_string": r.quality_string or "",
                    }
                    for r in fast5_reads
                ]

            if suffix in (".bam", ".cram"):
                from ..io.bam import read_long_read_bam

                bam_alignments = read_long_read_bam(reads_path)
                return [
                    {
                        "read_id": a.read_name,
                        "read_name": a.read_name,
                        "sequence": a.query_sequence,
                        "quality_string": "",
                        "reference_name": a.reference_name,
                        "reference_start": a.reference_start,
                        "reference_end": a.reference_end,
                        "mapping_quality": a.mapping_quality,
                        "cigar_string": a.cigar_string,
                        "cigar_tuples": a.cigar_tuples,
                        "is_supplementary": a.is_supplementary,
                        "is_secondary": a.is_secondary,
                        "is_reverse": a.is_reverse,
                        "is_unmapped": a.is_unmapped,
                        "tags": a.tags,
                        "query_length": a.query_length,
                    }
                    for a in bam_alignments
                ]

            raise ValueError(f"Unsupported input file format: {suffix}")

        raise TypeError(f"Unsupported input type: {type(reads).__name__}")

    def _load_alignments(self, context: dict[str, Any]) -> list[dict[str, Any]]:
        """Load alignments from step context.

        Args:
            context: Step context containing 'alignments' key.

        Returns:
            List of alignment dictionaries.
        """
        alignments = context.get("alignments", context)
        return self._load_alignments_sync(alignments)

    def _load_alignments_sync(self, alignments: Any) -> list[dict[str, Any]]:
        """Load alignments from various input sources.

        Args:
            alignments: List of dicts, LongReadAlignment objects, or Path to BAM.

        Returns:
            List of alignment dictionaries.
        """
        if isinstance(alignments, (list, tuple)):
            result: list[dict[str, Any]] = []
            for a in alignments:
                if isinstance(a, dict):
                    result.append(a)
                elif hasattr(a, "read_name"):
                    d: dict[str, Any] = {
                        "read_name": a.read_name,
                        "query_sequence": getattr(a, "query_sequence", ""),
                        "reference_name": getattr(a, "reference_name", ""),
                        "reference_start": getattr(a, "reference_start", 0),
                        "reference_end": getattr(a, "reference_end", 0),
                        "mapping_quality": getattr(a, "mapping_quality", 0),
                        "cigar_string": getattr(a, "cigar_string", ""),
                        "cigar_tuples": getattr(a, "cigar_tuples", []),
                        "is_supplementary": getattr(a, "is_supplementary", False),
                        "is_secondary": getattr(a, "is_secondary", False),
                        "is_reverse": getattr(a, "is_reverse", False),
                        "is_unmapped": getattr(a, "is_unmapped", False),
                        "tags": getattr(a, "tags", {}),
                        "query_length": getattr(a, "query_length", 0),
                    }
                    result.append(d)
            return result

        if isinstance(alignments, Path) or (isinstance(alignments, str) and Path(alignments).exists()):
            from ..io.bam import read_long_read_bam

            bam_alignments = read_long_read_bam(Path(alignments))
            return [
                {
                    "read_name": a.read_name,
                    "query_sequence": a.query_sequence,
                    "reference_name": a.reference_name,
                    "reference_start": a.reference_start,
                    "reference_end": a.reference_end,
                    "mapping_quality": a.mapping_quality,
                    "cigar_string": a.cigar_string,
                    "cigar_tuples": a.cigar_tuples,
                    "is_supplementary": a.is_supplementary,
                    "is_secondary": a.is_secondary,
                    "is_reverse": a.is_reverse,
                    "is_unmapped": a.is_unmapped,
                    "tags": a.tags,
                    "query_length": a.query_length,
                }
                for a in bam_alignments
            ]

        raise TypeError(f"Unsupported alignment input type: {type(alignments).__name__}")


# --- Module-level helpers ---


def _extract_methylation_features(
    reads: list[dict[str, Any]],
    mod_type: str = "5mC",
) -> list[dict[str, Any]]:
    """Extract methylation signal features from read dictionaries.

    Looks for pre-computed methylation features in the read metadata,
    or constructs synthetic features from alignment methylation tags.

    Args:
        reads: List of read dictionaries.
        mod_type: Modification type to extract ('5mC' or '6mA').

    Returns:
        List of feature dictionaries suitable for call_5mc() or call_6ma().
    """
    features: list[dict[str, Any]] = []

    for read in reads:
        # Check for pre-computed features
        read_features = read.get("methylation_features", [])
        if read_features:
            for feat in read_features:
                if feat.get("modification_type", mod_type) == mod_type:
                    features.append(feat)
            continue

        # Try to extract from methylation tags (BAM MM/ML tags)
        tags = read.get("tags", {})
        if isinstance(tags, dict):
            mm_tag = tags.get("MM") or tags.get("Mm")
            if mm_tag:
                from ..io.bam import _parse_methylation_from_tags

                parsed = _parse_methylation_from_tags(tags)
                modifications = parsed.get("modifications", [])
                for mod in modifications:
                    base = mod.get("base", "")
                    mod_code = mod.get("modification", "")

                    # Map modification codes to types
                    is_target = False
                    if mod_type == "5mC" and base == "C" and mod_code in ("m", "5mC"):
                        is_target = True
                    elif mod_type == "6mA" and base == "A" and mod_code in ("a", "6mA"):
                        is_target = True

                    if is_target:
                        deltas = mod.get("delta_positions", [])
                        probs = mod.get("probabilities", [])
                        chrom = read.get("reference_name", "")
                        ref_start = read.get("reference_start", 0)

                        for idx in range(len(deltas)):
                            position = ref_start + sum(deltas[: idx + 1]) + idx
                            prob = probs[idx] / 255.0 if idx < len(probs) else 0.0
                            features.append(
                                {
                                    "position": position,
                                    "chromosome": chrom,
                                    "strand": "+",
                                    "context_shift": prob * 5.0,  # Scale probability to approximate pA shift
                                    "dwell_time": 1.0,
                                    "current_std": 0.5,
                                    "context": "CpG" if mod_type == "5mC" else "A",
                                    "modification_type": mod_type,
                                }
                            )

    return features


def _write_json_safe(data: Any, path: Path) -> None:
    """Write data to JSON file, handling non-serializable types.

    Args:
        data: Data to serialize.
        path: Output path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    def _default(obj: Any) -> Any:
        if isinstance(obj, Path):
            return str(obj)
        if hasattr(obj, "__dataclass_fields__"):
            return {k: getattr(obj, k) for k in obj.__dataclass_fields__}
        if hasattr(obj, "__dict__"):
            return {k: v for k, v in obj.__dict__.items() if not k.startswith("_")}
        return str(obj)

    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=_default)
