"""Core orchestration engine for long-read analysis workflows.

Provides the LongReadOrchestrator class with core methods for pipeline execution,
dependency resolution, and data loading. Pipeline-specific methods (QC, assembly,
methylation, SV) are in pipeline_stages.py.

Each pipeline is defined as a directed acyclic graph (DAG) of PipelineStep
objects. The orchestrator performs topological sorting to determine execution
order and propagates results through a shared context dictionary.
"""

from __future__ import annotations

import json
import time
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

        Args:
            reads: List of read dictionaries or a Path to a FAST5/BAM file.

        Returns:
            PipelineResult with QC metrics and report paths.
        """
        from .pipeline_stages import _build_qc_steps

        steps = _build_qc_steps(self, reads)
        return self._run_steps("qc", steps)

    def run_assembly_pipeline(self, reads: list[dict[str, Any]] | Path) -> PipelineResult:
        """Run the assembly pipeline: QC, overlap, consensus, polish, stats.

        Args:
            reads: List of read dictionaries or Path to a FAST5/BAM file.

        Returns:
            PipelineResult with assembly contigs and statistics.
        """
        from .pipeline_stages import _build_assembly_steps

        steps = _build_assembly_steps(self, reads)
        return self._run_steps("assembly", steps)

    def run_methylation_pipeline(self, reads: list[dict[str, Any]] | Path) -> PipelineResult:
        """Run the methylation pipeline: load, call mods, aggregate, diff, plots.

        Args:
            reads: List of read dictionaries or a Path to a BAM file.

        Returns:
            PipelineResult with methylation calls and analysis results.
        """
        from .pipeline_stages import _build_methylation_steps

        steps = _build_methylation_steps(self, reads)
        return self._run_steps("methylation", steps)

    def run_sv_pipeline(self, alignments: list[dict[str, Any]] | Path) -> PipelineResult:
        """Run the SV pipeline: load, detect SVs, phase, annotate, plots.

        Args:
            alignments: List of alignment dictionaries or Path to a BAM file.

        Returns:
            PipelineResult with structural variant calls and statistics.
        """
        from .pipeline_stages import _build_sv_steps

        steps = _build_sv_steps(self, alignments)
        return self._run_steps("sv", steps)

    def run_full_pipeline(self, reads: list[dict[str, Any]] | Path) -> PipelineResult:
        """Run the full pipeline: QC, Assembly, Methylation, SV calling.

        Args:
            reads: List of read dictionaries or Path to input file.

        Returns:
            PipelineResult with combined results from all sub-pipelines.
        """
        from .pipeline_stages import _run_full_pipeline

        return _run_full_pipeline(self, reads)

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
