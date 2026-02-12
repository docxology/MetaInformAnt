"""Phenotype analysis pipeline orchestration.

Provides configurable multi-step pipelines for processing phenotype data
across behavioral, chemical, morphological, electronic, and sonic domains.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union

from metainformant.core.utils import logging
from metainformant.core.utils.errors import ValidationError

logger = logging.get_logger(__name__)


@dataclass
class PipelineConfig:
    """Configuration for a phenotype analysis pipeline.

    Attributes:
        name: Pipeline name for identification.
        phenotype_types: List of phenotype domains to analyze.
        input_path: Path to input data.
        output_path: Path for results.
        parameters: Domain-specific analysis parameters.
        steps: Ordered list of processing steps to execute.
    """

    name: str = "phenotype_pipeline"
    phenotype_types: List[str] = field(default_factory=lambda: ["morphological"])
    input_path: Optional[str] = None
    output_path: Optional[str] = None
    parameters: Dict[str, Any] = field(default_factory=dict)
    steps: List[str] = field(default_factory=lambda: ["load", "validate", "analyze", "summarize"])

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> PipelineConfig:
        """Load configuration from YAML file."""
        try:
            import yaml
        except ImportError:
            raise ImportError("YAML loading requires pyyaml: uv pip install pyyaml")

        with open(path) as f:
            data = yaml.safe_load(f)

        return cls(
            name=data.get("name", "phenotype_pipeline"),
            phenotype_types=data.get("phenotype_types", ["morphological"]),
            input_path=data.get("input_path"),
            output_path=data.get("output_path"),
            parameters=data.get("parameters", {}),
            steps=data.get("steps", ["load", "validate", "analyze", "summarize"]),
        )

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> PipelineConfig:
        """Load configuration from dictionary."""
        return cls(
            name=data.get("name", "phenotype_pipeline"),
            phenotype_types=data.get("phenotype_types", ["morphological"]),
            input_path=data.get("input_path"),
            output_path=data.get("output_path"),
            parameters=data.get("parameters", {}),
            steps=data.get("steps", ["load", "validate", "analyze", "summarize"]),
        )

    def validate(self) -> List[str]:
        """Validate configuration. Returns list of warnings/errors."""
        issues = []
        valid_types = {"morphological", "behavioral", "chemical", "electronic", "sonic"}
        for pt in self.phenotype_types:
            if pt not in valid_types:
                issues.append(f"Unknown phenotype type: {pt}")

        valid_steps = {"load", "validate", "preprocess", "analyze", "summarize", "export"}
        for step in self.steps:
            if step not in valid_steps:
                issues.append(f"Unknown step: {step}")

        return issues


@dataclass
class PipelineResult:
    """Result of a pipeline execution.

    Attributes:
        success: Whether pipeline completed successfully.
        config: The pipeline configuration used.
        outputs: Dict mapping step names to outputs.
        errors: List of errors encountered.
        warnings: List of warnings.
        metrics: Pipeline performance metrics.
        timestamp: When pipeline completed.
    """

    success: bool
    config: PipelineConfig
    outputs: Dict[str, Any] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    metrics: Dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary for serialization."""
        return {
            "success": self.success,
            "config_name": self.config.name,
            "phenotype_types": self.config.phenotype_types,
            "steps_executed": list(self.outputs.keys()),
            "errors": self.errors,
            "warnings": self.warnings,
            "metrics": self.metrics,
            "timestamp": self.timestamp,
        }

    def save_json(self, path: Union[str, Path]) -> None:
        """Save result summary to JSON file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2, default=str)


class PhenotypePipeline:
    """Orchestrates multi-step phenotype analysis pipelines.

    Supports configurable step sequences, domain-specific analyzers,
    and consistent error handling.

    Example:
        ```python
        config = PipelineConfig(
            name="ant_morphometry",
            phenotype_types=["morphological"],
            steps=["load", "validate", "analyze", "summarize"]
        )
        pipeline = PhenotypePipeline(config)
        result = pipeline.run(data)
        ```
    """

    def __init__(self, config: PipelineConfig) -> None:
        self.config = config
        self._step_functions: Dict[str, Callable] = {
            "load": self._step_load,
            "validate": self._step_validate,
            "preprocess": self._step_preprocess,
            "analyze": self._step_analyze,
            "summarize": self._step_summarize,
            "export": self._step_export,
        }
        self._state: Dict[str, Any] = {}

    def register_step(self, name: str, func: Callable) -> None:
        """Register a custom step function."""
        self._step_functions[name] = func

    def run(self, data: Optional[Any] = None) -> PipelineResult:
        """Execute the pipeline.

        Args:
            data: Optional input data. If None, attempts to load from config.input_path.

        Returns:
            PipelineResult with outputs and status.
        """
        logger.info(f"Starting pipeline: {self.config.name}")

        # Validate config (accept registered custom steps too)
        config_issues = self.config.validate()
        registered_steps = set(self._step_functions.keys())
        warnings = [w for w in config_issues if "warning" in w.lower()]
        errors = [
            e
            for e in config_issues
            if "warning" not in e.lower()
            and not (e.startswith("Unknown step:") and e.split(": ", 1)[1] in registered_steps)
        ]

        if errors:
            return PipelineResult(
                success=False,
                config=self.config,
                errors=errors,
                warnings=warnings,
            )

        # Initialize state
        self._state = {
            "raw_data": data,
            "processed_data": None,
            "results": {},
            "summary": {},
        }

        outputs: Dict[str, Any] = {}
        all_errors: List[str] = list(errors)
        all_warnings: List[str] = list(warnings)

        # Execute steps
        start_time = datetime.now()

        for step_name in self.config.steps:
            if step_name not in self._step_functions:
                all_errors.append(f"Step not found: {step_name}")
                continue

            logger.info(f"Executing step: {step_name}")
            try:
                step_result = self._step_functions[step_name]()
                outputs[step_name] = step_result

                # Collect step warnings
                if isinstance(step_result, dict) and "warnings" in step_result:
                    all_warnings.extend(step_result["warnings"])

            except Exception as e:
                error_msg = f"Step '{step_name}' failed: {str(e)}"
                logger.error(error_msg)
                all_errors.append(error_msg)
                break

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        metrics = {
            "duration_seconds": duration,
            "steps_completed": len(outputs),
            "steps_total": len(self.config.steps),
            "phenotype_types": self.config.phenotype_types,
        }

        success = len(all_errors) == 0 and len(outputs) == len(self.config.steps)

        result = PipelineResult(
            success=success,
            config=self.config,
            outputs=outputs,
            errors=all_errors,
            warnings=all_warnings,
            metrics=metrics,
        )

        logger.info(f"Pipeline completed: success={success}, duration={duration:.2f}s")
        return result

    def _step_load(self) -> Dict[str, Any]:
        """Load input data."""
        if self._state["raw_data"] is not None:
            return {"status": "data_provided", "n_items": self._count_items(self._state["raw_data"])}

        if self.config.input_path:
            path = Path(self.config.input_path)
            if path.suffix == ".json":
                with open(path) as f:
                    self._state["raw_data"] = json.load(f)
            else:
                raise ValidationError(f"Unsupported input format: {path.suffix}")

            return {"status": "loaded", "path": str(path), "n_items": self._count_items(self._state["raw_data"])}

        return {"status": "no_data", "warnings": ["No input data provided"]}

    def _step_validate(self) -> Dict[str, Any]:
        """Validate input data."""
        data = self._state["raw_data"]
        if data is None:
            return {"status": "skipped", "reason": "no_data"}

        issues = []
        validated_count = 0

        # Type-specific validation
        if isinstance(data, list):
            for i, item in enumerate(data):
                if not isinstance(item, dict):
                    issues.append(f"Item {i} is not a dict")
                else:
                    validated_count += 1
        elif isinstance(data, dict):
            validated_count = 1
        else:
            issues.append(f"Data must be dict or list, got {type(data).__name__}")

        self._state["processed_data"] = data
        return {
            "status": "validated" if not issues else "warnings",
            "validated_count": validated_count,
            "warnings": issues,
        }

    def _step_preprocess(self) -> Dict[str, Any]:
        """Preprocess data for analysis."""
        data = self._state.get("processed_data") or self._state.get("raw_data")
        if data is None:
            return {"status": "skipped", "reason": "no_data"}

        # Apply preprocessing based on phenotype type
        params = self.config.parameters.get("preprocess", {})
        normalize = params.get("normalize", False)
        filter_missing = params.get("filter_missing", True)

        processed = data
        ops_applied = []

        if filter_missing and isinstance(processed, list):
            original_len = len(processed)
            processed = [d for d in processed if d]
            if len(processed) < original_len:
                ops_applied.append(f"filtered_empty:{original_len - len(processed)}")

        self._state["processed_data"] = processed
        return {"status": "preprocessed", "operations": ops_applied}

    def _step_analyze(self) -> Dict[str, Any]:
        """Run domain-specific analysis."""
        data = self._state.get("processed_data") or self._state.get("raw_data")
        if data is None:
            return {"status": "skipped", "reason": "no_data"}

        results = {}
        for ptype in self.config.phenotype_types:
            analyzer = self._get_analyzer(ptype)
            if analyzer:
                try:
                    result = analyzer(data, self.config.parameters.get(ptype, {}))
                    results[ptype] = result
                except Exception as e:
                    results[ptype] = {"error": str(e)}

        self._state["results"] = results
        return {"status": "analyzed", "domains": list(results.keys()), "results": results}

    def _step_summarize(self) -> Dict[str, Any]:
        """Generate summary statistics."""
        results = self._state.get("results", {})

        summary = {
            "n_domains": len(results),
            "domains": list(results.keys()),
            "per_domain": {},
        }

        for domain, data in results.items():
            if isinstance(data, dict):
                summary["per_domain"][domain] = {
                    "n_keys": len(data),
                    "has_error": "error" in data,
                }

        self._state["summary"] = summary
        return {"status": "summarized", "summary": summary}

    def _step_export(self) -> Dict[str, Any]:
        """Export results to output path."""
        if not self.config.output_path:
            return {"status": "skipped", "reason": "no_output_path"}

        output_path = Path(self.config.output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        export_data = {
            "config": self.config.name,
            "results": self._state.get("results", {}),
            "summary": self._state.get("summary", {}),
            "timestamp": datetime.now().isoformat(),
        }

        with open(output_path, "w") as f:
            json.dump(export_data, f, indent=2, default=str)

        return {"status": "exported", "path": str(output_path)}

    def _get_analyzer(self, phenotype_type: str) -> Optional[Callable]:
        """Get the analyzer function for a phenotype type."""
        analyzers = {
            "morphological": self._analyze_morphological,
            "behavioral": self._analyze_behavioral,
            "chemical": self._analyze_chemical,
            "electronic": self._analyze_electronic,
            "sonic": self._analyze_sonic,
        }
        return analyzers.get(phenotype_type)

    def _analyze_morphological(self, data: Any, params: Dict) -> Dict[str, Any]:
        """Morphological phenotype analysis."""
        from ..morphological.measurement import Measurement
        from ..morphological.profile import MorphometricProfile, summary_statistics

        if isinstance(data, list) and all(isinstance(d, MorphometricProfile) for d in data):
            return summary_statistics(data)

        # Basic analysis for raw dict data
        if isinstance(data, list):
            return {"n_specimens": len(data), "type": "raw_data"}
        return {"type": "single_specimen"}

    def _analyze_behavioral(self, data: Any, params: Dict) -> Dict[str, Any]:
        """Behavioral phenotype analysis."""
        from ..behavior.sequence import BehaviorSequence

        if isinstance(data, list) and all(isinstance(d, BehaviorSequence) for d in data):
            diversities = [seq.shannon_diversity() for seq in data]
            return {
                "n_sequences": len(data),
                "mean_diversity": sum(diversities) / len(diversities) if diversities else 0,
            }
        return {"n_items": len(data) if isinstance(data, list) else 1}

    def _analyze_chemical(self, data: Any, params: Dict) -> Dict[str, Any]:
        """Chemical phenotype analysis."""
        from ..chemical.profile import ChemicalProfile

        if isinstance(data, list) and all(isinstance(d, ChemicalProfile) for d in data):
            diversities = [p.shannon_diversity() for p in data]
            return {
                "n_profiles": len(data),
                "mean_diversity": sum(diversities) / len(diversities) if diversities else 0,
            }
        return {"n_items": len(data) if isinstance(data, list) else 1}

    def _analyze_electronic(self, data: Any, params: Dict) -> Dict[str, Any]:
        """Electronic/tracking phenotype analysis."""
        from ..electronic.tracking import Trajectory

        if isinstance(data, list) and all(isinstance(d, Trajectory) for d in data):
            distances = [t.total_distance() for t in data]
            return {
                "n_trajectories": len(data),
                "mean_distance": sum(distances) / len(distances) if distances else 0,
            }
        return {"n_items": len(data) if isinstance(data, list) else 1}

    def _analyze_sonic(self, data: Any, params: Dict) -> Dict[str, Any]:
        """Sonic phenotype analysis."""
        from ..sonic.signal import AcousticSignal

        if isinstance(data, list) and all(isinstance(d, AcousticSignal) for d in data):
            frequencies = [s.dominant_frequency() for s in data]
            return {
                "n_signals": len(data),
                "mean_dominant_freq": sum(frequencies) / len(frequencies) if frequencies else 0,
            }
        return {"n_items": len(data) if isinstance(data, list) else 1}

    @staticmethod
    def _count_items(data: Any) -> int:
        """Count items in data."""
        if isinstance(data, list):
            return len(data)
        elif isinstance(data, dict):
            return len(data)
        return 1
