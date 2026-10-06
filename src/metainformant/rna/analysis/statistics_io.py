"""Strict JSON deserialization for predeclared RNA analysis contracts."""

from __future__ import annotations
import json
from pathlib import Path
from metainformant.rna.analysis.statistics_contract import (
    AnalysisProvenance,
    PredeclaredDesign,
    SensitivityAnalysis,
)


def read_analysis_provenance(path: Path) -> AnalysisProvenance:
    """Deserialize an ``AnalysisProvenance`` contract from a JSON file.

    Fail-closed: unknown fields, sensitivity entries that are not objects,
    and constructor mismatches are all rejected before any gate runs.

    Raises:
        ValueError: If the file is unreadable, is not a JSON object,
            declares unknown fields, or does not match the record shape.
    """

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"contract file {path} is unreadable: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"contract file {path} must contain a JSON object")
    known = set(AnalysisProvenance.__dataclass_fields__)
    unknown = sorted(set(payload) - known)
    if unknown:
        raise ValueError(f"contract file {path} declares unknown fields: {unknown}")
    kwargs = dict(payload)
    entries_payload = kwargs.pop("sensitivity_analyses", None)
    if entries_payload is None:
        entries_payload = []
    if not isinstance(entries_payload, list):
        raise ValueError(
            f"contract file {path} declares sensitivity_analyses as "
            f"{type(entries_payload).__name__}; a JSON list is required"
        )
    entries: list[SensitivityAnalysis] = []
    for entry in entries_payload:
        if not isinstance(entry, dict):
            raise ValueError(
                f"contract file {path} declares a non-object sensitivity entry: {entry!r}"
            )
        try:
            entry_record = SensitivityAnalysis(
                **{**entry, "varied_values": tuple(entry.get("varied_values", ()))}
            )
        except TypeError as exc:
            raise ValueError(
                f"contract file {path} declares a sensitivity entry that does not "
                f"match the SensitivityAnalysis record: {exc}"
            ) from exc
        entries.append(entry_record)
    kwargs["sensitivity_analyses"] = tuple(entries)
    design_payload = kwargs.pop("design_declaration", None)
    if design_payload is not None:
        if not isinstance(design_payload, dict):
            raise ValueError(
                f"contract file {path} declares design_declaration as "
                f"{type(design_payload).__name__}; a JSON object is required"
            )
        unknown_design_fields = sorted(set(design_payload) - {"covariate_strata"})
        if unknown_design_fields:
            raise ValueError(
                f"contract file {path} declares unknown design_declaration fields: "
                f"{unknown_design_fields}"
            )
        strata_payload = design_payload.get("covariate_strata")
        if not isinstance(strata_payload, dict) or not strata_payload:
            raise ValueError(
                f"contract file {path} declares design_declaration.covariate_strata as "
                f"{strata_payload!r}; a non-empty JSON object is required"
            )
        covariate_strata: dict[str, tuple[str, ...]] = {}
        for covariate, levels in strata_payload.items():
            if not isinstance(levels, list) or not levels:
                raise ValueError(
                    f"contract file {path} declares strata for covariate {covariate!r} as "
                    f"{levels!r}; a non-empty JSON list is required"
                )
            if not all(isinstance(level, str) for level in levels):
                raise ValueError(
                    f"contract file {path} declares a non-string stratum level for "
                    f"covariate {covariate!r}: {levels!r}"
                )
            covariate_strata[str(covariate)] = tuple(levels)
        kwargs["design_declaration"] = PredeclaredDesign(
            covariate_strata=covariate_strata
        )
    try:
        contract = AnalysisProvenance(**kwargs)
    except TypeError as exc:
        raise ValueError(
            f"contract file {path} does not match the AnalysisProvenance record: {exc}"
        ) from exc
    return contract
