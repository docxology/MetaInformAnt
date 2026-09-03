"""Predeclared statistical contract for cross-species RNA analyses.

Implements the analysis-provenance and descriptive/inferential boundary
required by
``projects/hymenoptera_amalgkit/docs/manuscript/statistical_analysis_plan.md``
(sections 1, 4, 6, 8, and 9):

- :class:`AnalysisProvenance` is the structured record declared before an
  analysis runs: analysis identifier, estimand, biological replicate unit,
  random seed, resampling count, null model, multiple-testing procedure
  (family, method, tested-feature count), analysis role, and software
  versions.
- Descriptive outputs (the fingerprint divergence matrix and the
  feature-resampling sensitivity table in
  :mod:`metainformant.rna.analysis.cross_species`) carry
  ``attrs["role"] == "descriptive"``. They are permutation/sensitivity
  scores, never p-values, and never confidence intervals.
- The only inferential path is :func:`declared_inferential_bh_fdr`, which is
  GATED (mirroring ``wilcoxon_duplication_specificity`` in
  ``tissue_specificity.py``): it refuses to run unless the caller passes an
  already-validated contract declared ``inferential`` AND sets
  ``evidence_manifest_frozen=True``. It applies the declared
  Benjamini-Hochberg FDR procedure before returning adjusted values.
- :func:`validate_orthology_profile_invariants` and
  :func:`validate_species_tree_invariants` fail closed with explicit error
  types on orthology-bridge and species-tree violations (plan sections 4
  and 8).

Every validator fails closed: a record with missing or placeholder fields,
or an invariant violation, raises before any result can be produced.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import pandas as pd

__all__ = [
    "AnalysisProvenance",
    "DESCRIPTIVE_ROLE",
    "INFERENTIAL_ROLE",
    "OrthologyInvariantError",
    "ProvenanceError",
    "StatisticsContractError",
    "TreeInvariantError",
    "benjamini_hochberg_fdr",
    "declared_inferential_bh_fdr",
    "render_analysis_provenance_block",
    "result_role",
    "validate_analysis_provenance",
    "validate_orthology_profile_invariants",
    "validate_species_tree_invariants",
]

DESCRIPTIVE_ROLE = "descriptive"
INFERENTIAL_ROLE = "inferential"

_ALLOWED_ROLES = frozenset({DESCRIPTIVE_ROLE, INFERENTIAL_ROLE})
# Predeclared multiplicity procedures (statistical_analysis_plan.md section 6).
_ALLOWED_MT_METHODS = frozenset({"bh-fdr", "benjamini-hochberg", "bonferroni"})
# Strings that make a declared field a placeholder rather than a declaration.
_PLACEHOLDER_STRINGS = frozenset(
    {"", "na", "n/a", "none", "null", "todo", "tbd", "placeholder", "unknown", "pending", "?"}
)


class StatisticsContractError(ValueError):
    """Base class for fail-closed statistical-contract violations."""


class ProvenanceError(StatisticsContractError):
    """An analysis-provenance record is missing or contains placeholder fields."""


class OrthologyInvariantError(StatisticsContractError):
    """An orthology x species presence table violates a declared invariant."""


class TreeInvariantError(StatisticsContractError):
    """A species tree is unrooted, malformed, or has conflicting labels."""


@dataclass(frozen=True)
class AnalysisProvenance:
    """Predeclared analysis-provenance record (plan sections 1, 4, 6, 9).

    The record is frozen: after declaration it is the immutable contract the
    analysis runs under. Descriptive analyses declare
    ``analysis_role="descriptive"``; only contracts declared
    ``analysis_role="inferential"`` with a BH-FDR procedure may be used with
    :func:`declared_inferential_bh_fdr`.
    """

    analysis_id: str
    estimand: str
    replicate_unit: str
    random_seed: int
    resampling_count: int
    null_model: str
    multiple_testing_family: str | None
    multiple_testing_method: str | None
    tested_feature_count: int | None
    software_versions: Mapping[str, str]
    analysis_role: str = DESCRIPTIVE_ROLE


def _require_declared(value: object, field: str) -> str:
    """Reject missing or placeholder strings for a declared field."""
    if not isinstance(value, str) or value.strip().lower() in _PLACEHOLDER_STRINGS:
        raise ProvenanceError(
            f"analysis provenance field '{field}' is missing or a placeholder "
            f"(got {value!r}); a declared value is required"
        )
    return value


def validate_analysis_provenance(record: AnalysisProvenance) -> None:
    """Fail closed on missing, placeholder, or inconsistent provenance fields.

    Raises:
        ProvenanceError: On any missing/placeholder field, a negative seed,
            a role mismatch, or an empty software-version mapping. The
            multiplicity fields are role-conditional and fail closed in both
            directions: ``multiple_testing_family``,
            ``multiple_testing_method``, and ``tested_feature_count`` must
            each be None (or the literal ``'not-applicable'`` for the two
            string fields) exactly when ``analysis_role`` is
            ``'descriptive'`` (no inferential test was performed); a real
            declared family, BH-FDR | Benjamini-Hochberg | Bonferroni
            method, and positive tested-feature count are required exactly
            when the role is ``'inferential'``.
        TypeError: If ``record`` is not an :class:`AnalysisProvenance`.
    """
    if not isinstance(record, AnalysisProvenance):
        raise TypeError(
            "validate_analysis_provenance requires an AnalysisProvenance record; "
            "ad-hoc dictionaries cannot bind an analysis to its predeclared contract"
        )
    for field in (
        "analysis_id",
        "estimand",
        "replicate_unit",
        "null_model",
    ):
        _require_declared(getattr(record, field), field)

    if record.analysis_role not in _ALLOWED_ROLES:
        raise ProvenanceError(f"analysis_role must be one of {sorted(_ALLOWED_ROLES)}, got {record.analysis_role!r}")

    if not isinstance(record.random_seed, int) or isinstance(record.random_seed, bool) or record.random_seed < 0:
        raise ProvenanceError(f"random_seed must be a non-negative integer, got {record.random_seed!r}")
    if (
        not isinstance(record.resampling_count, int)
        or isinstance(record.resampling_count, bool)
        or record.resampling_count < 1
    ):
        raise ProvenanceError(f"resampling_count must be a positive integer, got {record.resampling_count!r}")
    if record.tested_feature_count is None:
        if record.analysis_role != DESCRIPTIVE_ROLE:
            raise ProvenanceError(
                "analysis_role='inferential' must declare the number of tested "
                "features; got tested_feature_count=None"
            )
    elif record.analysis_role == DESCRIPTIVE_ROLE:
        raise ProvenanceError(
            f"analysis_role='descriptive' tests no features, so declaring "
            f"tested_feature_count={record.tested_feature_count!r} misstates the "
            f"analysis; use None (rendered 'not-applicable')"
        )
    elif (
        not isinstance(record.tested_feature_count, int)
        or isinstance(record.tested_feature_count, bool)
        or record.tested_feature_count < 1
    ):
        raise ProvenanceError(
            f"tested_feature_count must be a positive integer for inferential "
            f"analyses, got {record.tested_feature_count!r}"
        )

    # Multiplicity procedures belong to inferential tests only (plan section
    # 6): a descriptive lane applies no procedure, so recording one would
    # misstate the analysis. The symmetry is enforced fail-closed in BOTH
    # directions.
    method = record.multiple_testing_method
    declared_not_applicable = method is None or (isinstance(method, str) and method.strip().lower() == "not-applicable")
    if declared_not_applicable:
        if record.analysis_role != DESCRIPTIVE_ROLE:
            raise ProvenanceError(
                f"analysis_role='inferential' must declare a real multiple-testing "
                f"procedure; got multiple_testing_method={method!r}"
            )
    else:
        if record.analysis_role != INFERENTIAL_ROLE:
            raise ProvenanceError(
                f"analysis_role='descriptive' applies no multiple-testing procedure, so "
                f"declaring multiple_testing_method={method!r} misstates the analysis; use "
                f"None or 'not-applicable' for descriptive lanes"
            )
        declared = _require_declared(method, "multiple_testing_method").strip().lower()
        if declared not in _ALLOWED_MT_METHODS:
            raise ProvenanceError(
                f"multiple_testing_method must be one of {sorted(_ALLOWED_MT_METHODS)}, got "
                f"{record.multiple_testing_method!r}"
            )

    # The multiplicity FAMILY names the tested units of an inferential test
    # (plan section 6); a descriptive lane has none, so declaring one would
    # misstate the analysis. Same fail-closed symmetry as the method.
    family = record.multiple_testing_family
    family_not_applicable = family is None or (isinstance(family, str) and family.strip().lower() == "not-applicable")
    if family_not_applicable:
        if record.analysis_role != DESCRIPTIVE_ROLE:
            raise ProvenanceError(
                f"analysis_role='inferential' must declare a real multiple-testing "
                f"family; got multiple_testing_family={family!r}"
            )
    else:
        if record.analysis_role != INFERENTIAL_ROLE:
            raise ProvenanceError(
                f"analysis_role='descriptive' performs no inferential test, so "
                f"declaring multiple_testing_family={family!r} misstates the analysis; "
                f"use None or 'not-applicable' for descriptive lanes"
            )
        _require_declared(family, "multiple_testing_family")

    if not isinstance(record.software_versions, Mapping) or not record.software_versions:
        raise ProvenanceError("software_versions must be a non-empty mapping of package name to version")
    for name, version in record.software_versions.items():
        _require_declared(name, "software_versions key")
        _require_declared(version, f"software_versions[{name!r}]")


def render_analysis_provenance_block(record: AnalysisProvenance) -> list[str]:
    """Render the provenance record as additive ``key: value`` summary lines.

    The lines use the same ``key: value`` shape as the project's
    ``analysis_summary.txt`` (written by ``run_cross_species_analysis.py``),
    so they can be appended without changing existing keys. Validation runs
    first: a record that would render placeholder provenance fails closed.

    Returns:
        Deterministic list of ``analysis_provenance_*`` lines.
    """
    validate_analysis_provenance(record)
    versions = "; ".join(f"{name}={record.software_versions[name]}" for name in sorted(record.software_versions))
    return [
        f"analysis_provenance_role: {record.analysis_role}",
        f"analysis_provenance_analysis_id: {record.analysis_id}",
        f"analysis_provenance_estimand: {record.estimand}",
        f"analysis_provenance_replicate_unit: {record.replicate_unit}",
        f"analysis_provenance_random_seed: {record.random_seed}",
        f"analysis_provenance_resampling_count: {record.resampling_count}",
        f"analysis_provenance_null_model: {record.null_model}",
        "analysis_provenance_multiple_testing_family: "
        + ("not-applicable" if record.multiple_testing_family is None else record.multiple_testing_family),
        "analysis_provenance_multiple_testing_method: "
        + (
            "not-applicable"
            if record.multiple_testing_method is None
            else record.multiple_testing_method.strip().lower()
        ),
        "analysis_provenance_tested_feature_count: "
        + ("not-applicable" if record.tested_feature_count is None else str(record.tested_feature_count)),
        f"analysis_provenance_software_versions: {versions}",
    ]


def result_role(result: Any) -> str:
    """Return the declared role of a result object, failing closed.

    Descriptive results produced by ``cross_species.compute_fingerprint_*``
    carry ``attrs["role"] == "descriptive"``. Any result without an explicit
    role marker is refused rather than being silently treated as either
    descriptive or inferential.
    """
    attrs = getattr(result, "attrs", None)
    if not isinstance(attrs, Mapping) or "role" not in attrs:
        raise StatisticsContractError("result carries no declared role; refusing to classify unlabeled output")
    role = attrs["role"]
    if role not in _ALLOWED_ROLES:
        raise StatisticsContractError(f"result declares unknown role {role!r}")
    return str(role)


def benjamini_hochberg_fdr(p_values: Sequence[float]) -> list[float]:
    """Apply the Benjamini-Hochberg FDR procedure to raw p-values.

    Args:
        p_values: Non-empty sequence of raw p-values in [0, 1].

    Returns:
        Adjusted q-values in the same order as the input.

    Raises:
        StatisticsContractError: If the sequence is empty or any value is
            not a finite p-value in [0, 1].
    """
    values = list(p_values)
    if not values:
        raise StatisticsContractError("benjamini_hochberg_fdr requires at least one p-value")
    for index, value in enumerate(values):
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise StatisticsContractError(f"p-values must be numeric; index {index} got {value!r}")
        if value != value or value in (float("inf"), float("-inf")):
            raise StatisticsContractError(f"p-value at index {index} is not finite")
        if value < 0.0 or value > 1.0:
            raise StatisticsContractError(f"p-value at index {index} is outside [0, 1]: {value!r}")

    n = len(values)
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    adjusted = [0.0] * n
    running_min = 1.0
    for rank_index in range(n - 1, -1, -1):
        original_index, p_value = indexed[rank_index]
        rank = rank_index + 1
        running_min = min(running_min, p_value * n / rank)
        adjusted[original_index] = min(running_min, 1.0)
    return adjusted


def declared_inferential_bh_fdr(
    p_values: Sequence[float],
    contract: AnalysisProvenance,
    evidence_manifest_frozen: bool = False,
) -> dict[str, Any]:
    """GATED: apply the contract's declared BH-FDR procedure to raw p-values.

    POST-FREEZE USE ONLY. Descriptive permutation and sensitivity scores must
    never reach this function; inferential adjustment requires (a) an explicit
    ``evidence_manifest_frozen=True`` affirmation at every call site, and
    (b) a validated :class:`AnalysisProvenance` declared
    ``analysis_role="inferential"`` with a BH-FDR procedure whose
    ``tested_feature_count`` matches the supplied family exactly. Default
    ``evidence_manifest_frozen=False`` always refuses.

    Returns:
        Dict with ``role: "inferential"``, the declared method, raw and
        BH-adjusted p-values, and the contract fields echoed for provenance.

    Raises:
        RuntimeError: If ``evidence_manifest_frozen`` is not True.
        ProvenanceError: If the contract itself fails validation.
        StatisticsContractError: If the contract is not declared inferential
            with BH-FDR, or the tested-feature count does not match.
    """
    if not evidence_manifest_frozen:
        raise RuntimeError(
            "declared_inferential_bh_fdr is gated for post-freeze use: "
            "refusing to run while the evidence manifest is unfrozen"
        )
    validate_analysis_provenance(contract)
    if contract.analysis_role != INFERENTIAL_ROLE:
        raise StatisticsContractError(
            f"declared_inferential_bh_fdr requires a contract declared "
            f"analysis_role={INFERENTIAL_ROLE!r}, got {contract.analysis_role!r}"
        )
    if contract.multiple_testing_method.strip().lower() not in {"bh-fdr", "benjamini-hochberg"}:
        raise StatisticsContractError(
            f"declared_inferential_bh_fdr requires a BH-FDR procedure, got " f"{contract.multiple_testing_method!r}"
        )

    values = [float(value) for value in p_values]
    if contract.tested_feature_count != len(values):
        raise StatisticsContractError(
            f"contract declares tested_feature_count={contract.tested_feature_count} but "
            f"{len(values)} p-values were supplied; refusing to adjust a different family"
        )
    adjusted = benjamini_hochberg_fdr(values)
    return {
        "role": INFERENTIAL_ROLE,
        "analysis_id": contract.analysis_id,
        "estimand": contract.estimand,
        "replicate_unit": contract.replicate_unit,
        "random_seed": contract.random_seed,
        "null_model": contract.null_model,
        "multiple_testing_family": contract.multiple_testing_family,
        "multiple_testing_method": "bh-fdr",
        "tested_feature_count": len(values),
        "raw_p_values": values,
        "adjusted_p_values": adjusted,
        "gate": "post-freeze",
    }


# =============================================================================
# Orthology and species-tree invariants (plan sections 4 and 8)
# =============================================================================


def validate_orthology_profile_invariants(
    presence: pd.DataFrame,
    *,
    min_species_fraction: float = 0.5,
    min_species_per_orthogroup: int = 2,
) -> None:
    """Fail closed on orthology x species presence-table violations.

    The presence table has orthogroups as rows and species as columns; a
    truthy value means the orthogroup is present (mapped) in that species.
    Enforced invariants:

    - unique, non-placeholder orthogroup and species labels;
    - presence values restricted to boolean or numeric 0/1 — missing or
      ambiguous states must be explicit, never silently coerced (section 8);
    - every orthogroup maps to at least
      ``max(min_species_per_orthogroup, ceil(min_species_fraction * n_species))``
      species; violations are reported together with explicit reason codes.

    Raises:
        OrthologyInvariantError: On any violated invariant.
    """
    if not isinstance(presence, pd.DataFrame):
        raise OrthologyInvariantError("presence table must be a pandas DataFrame")
    if presence.empty:
        raise OrthologyInvariantError("presence table is empty")
    if presence.shape[1] < 2:
        raise OrthologyInvariantError(f"presence table must cover at least 2 species, got {presence.shape[1]}")
    if not 0.0 <= min_species_fraction <= 1.0:
        raise OrthologyInvariantError(f"min_species_fraction must be within [0, 1], got {min_species_fraction!r}")
    if (
        not isinstance(min_species_per_orthogroup, int)
        or isinstance(min_species_per_orthogroup, bool)
        or min_species_per_orthogroup < 1
    ):
        raise OrthologyInvariantError(
            f"min_species_per_orthogroup must be a positive integer, got " f"{min_species_per_orthogroup!r}"
        )

    duplicated_orthogroups = sorted({label for label in presence.index[presence.index.duplicated(keep=False)].tolist()})
    if duplicated_orthogroups:
        raise OrthologyInvariantError(f"duplicate orthogroup labels: {duplicated_orthogroups}")
    duplicated_species = sorted({label for label in presence.columns[presence.columns.duplicated(keep=False)].tolist()})
    if duplicated_species:
        raise OrthologyInvariantError(f"duplicate species labels: {duplicated_species}")

    bad_labels = [
        label
        for label in list(presence.index) + list(presence.columns)
        if not isinstance(label, str) or label.strip() == "" or label.strip().lower() in _PLACEHOLDER_STRINGS
    ]
    if bad_labels:
        raise OrthologyInvariantError(f"presence table contains missing or placeholder labels: {bad_labels}")

    non_numeric = [
        str(column)
        for column, dtype in presence.dtypes.items()
        if not pd.api.types.is_bool_dtype(dtype) and not pd.api.types.is_numeric_dtype(dtype)
    ]
    if non_numeric:
        raise OrthologyInvariantError(
            "presence values must be boolean or numeric 0/1; non-numeric columns: " + ", ".join(non_numeric)
        )
    numeric = presence.astype(float)
    if numeric.isna().any().any():
        missing = numeric.isna().any(axis=1)
        offenders = [str(label) for label in numeric.index[missing].tolist()]
        raise OrthologyInvariantError(
            "presence table contains missing values; absence must be encoded explicitly as 0, "
            "offending orthogroups: " + ", ".join(offenders)
        )
    unexpected_mask = (numeric != 0.0) & (numeric != 1.0)
    if unexpected_mask.any().any():
        offenders = [str(label) for label in numeric.index[unexpected_mask.any(axis=1)].tolist()]
        raise OrthologyInvariantError(
            "presence values outside {0, 1} are ambiguous and refused, offending orthogroups: " + ", ".join(offenders)
        )

    n_species = numeric.shape[1]
    required = max(min_species_per_orthogroup, math.ceil(min_species_fraction * n_species))
    counts = numeric.sum(axis=1)
    insufficient = counts[counts < required]
    if not insufficient.empty:
        detail = ", ".join(f"{label} ({int(count)}/{n_species} species)" for label, count in insufficient.items())
        raise OrthologyInvariantError(
            f"orthology bridge below declared coverage threshold "
            f"(required {required} of {n_species} species per orthogroup): {detail}"
        )


def _parse_newick(newick: str) -> dict[str, Any]:
    """Parse a plain Newick string (names, branch lengths) into nested dicts.

    Supports unquoted names, optional branch lengths, and optional internal
    labels. Any parse failure raises :class:`TreeInvariantError`.
    """
    text = newick.strip()
    if not text.endswith(";"):
        raise TreeInvariantError("Newick string must terminate with ';'")
    text = text[:-1]
    position = 0

    def parse_node() -> dict[str, Any]:
        nonlocal position
        children: list[dict[str, Any]] = []
        if position < len(text) and text[position] == "(":
            position += 1
            children.append(parse_node())
            while position < len(text) and text[position] == ",":
                position += 1
                children.append(parse_node())
            if position >= len(text) or text[position] != ")":
                raise TreeInvariantError(f"unbalanced parentheses at character {position}")
            position += 1
        start = position
        while position < len(text) and text[position] not in "(),:":
            position += 1
        name = text[start:position].strip()
        if position < len(text) and text[position] == ":":
            position += 1
            length_start = position
            while position < len(text) and text[position] not in "(),":
                position += 1
            length_text = text[length_start:position].strip()
            try:
                float(length_text)
            except ValueError:
                raise TreeInvariantError(f"invalid branch length {length_text!r} in Newick string") from None
        return {"name": name, "children": children}

    root = parse_node()
    if position != len(text):
        raise TreeInvariantError(f"unexpected characters after Newick root at position {position}")
    return root


def validate_species_tree_invariants(
    tree: str | Mapping[str, Any],
    *,
    rooted: bool | None = None,
    require_rooted: bool = True,
    require_bifurcating_root: bool = False,
) -> None:
    """Fail closed on malformed, unlabeled, or undeclared species trees.

    Accepts either a Newick string or the nested-dict representation used by
    ``cross_species.phylogenetic_expression_profile`` (keys ``name``,
    ``children``, ``distance``).

    Rootedness is CALLER-DECLARED PROVENANCE (statistical_analysis_plan.md
    section 5.3 requires the tree source and its provenance to be recorded).
    Plain Newick topology cannot establish biological rootedness: an
    unrooted tree can be serialized behind a bifurcating root, and a
    polytomous root is a valid hard polytomy in a rooted tree. Therefore the
    number of root children is never used to infer rootedness. When
    ``require_rooted`` is set, the caller must declare rootedness through
    ``rooted``:

    - ``rooted=None`` (default) fails closed with :class:`ProvenanceError`:
      an explicit rootedness declaration is required;
    - ``rooted=False`` raises :class:`TreeInvariantError`;
    - ``rooted=True`` passes.

    Enforced invariants:

    - every node is well-formed (leaves and, for dict input, internal nodes
      and the root carry non-placeholder names; internal names may be empty
      only in Newick input);
    - leaf labels are unique (duplicated taxa fail closed);
    - at least two leaves;
    - declared rootedness per ``require_rooted``;
    - ``require_bifurcating_root=True`` additionally asserts the structural
      property that the root node has exactly 2 children (a naming of what
      is actually checked; it does not prove biological rootedness).

    Raises:
        TreeInvariantError: On any violated structural invariant, malformed
            input, or a declared-unrooted tree when ``require_rooted`` is set.
        ProvenanceError: If ``require_rooted`` is set and ``rooted`` was not
            explicitly declared by the caller.
    """
    if isinstance(tree, str):
        root = _parse_newick(tree)
        internal_names_optional = True
    elif isinstance(tree, Mapping):
        root = tree
        internal_names_optional = False
    else:
        raise TreeInvariantError("species tree must be a Newick string or a nested-dict tree mapping")

    leaves: list[str] = []
    seen_leaves: dict[str, int] = {}
    root_children: list[Any] | None = None

    def visit(node: Any, depth: int, *, is_root: bool) -> None:
        nonlocal root_children
        if not isinstance(node, Mapping) or "name" not in node:
            raise TreeInvariantError("every tree node must be a mapping with a 'name' key")
        name = node["name"]
        if not isinstance(name, str):
            raise TreeInvariantError("tree node labels must be strings")
        label = name.strip()
        # Dict trees must name every node (phylogenetic_expression_profile
        # reads node["name"]); Newick legitimately allows unnamed internal
        # nodes and roots, so only leaves are held to the named standard.
        if not internal_names_optional and (label == "" or label.lower() in _PLACEHOLDER_STRINGS):
            raise TreeInvariantError(
                f"tree {'root' if is_root else 'internal'} label is missing or a placeholder: {name!r}"
            )
        children = node.get("children", [])
        if not isinstance(children, list):
            raise TreeInvariantError(f"node {name!r} has a non-list 'children' field")
        if is_root:
            root_children = children
        if not children:
            if not isinstance(name, str) or name.strip() == "" or name.strip().lower() in _PLACEHOLDER_STRINGS:
                raise TreeInvariantError(f"tree leaf label is missing or a placeholder: {name!r}")
            leaves.append(name)
            seen_leaves[name] = seen_leaves.get(name, 0) + 1
            return
        for child in children:
            visit(child, depth + 1, is_root=False)

    visit(root, 0, is_root=True)

    duplicates = sorted(label for label, count in seen_leaves.items() if count > 1)
    if duplicates:
        raise TreeInvariantError(f"species tree contains duplicate leaf labels: {duplicates}")
    if len(leaves) < 2:
        raise TreeInvariantError(f"species tree must have at least 2 leaves, got {len(leaves)}")
    assert root_children is not None
    if len(root_children) < 2:
        raise TreeInvariantError(f"species tree root has {len(root_children)} children; at least 2 are required")
    if require_rooted:
        # Rootedness cannot be inferred from topology; demand explicit
        # provenance from the caller (plan section 5.3).
        if rooted is None:
            raise ProvenanceError(
                "require_rooted is set but rootedness was not declared; pass rooted=True or "
                "rooted=False from recorded tree provenance (topology alone cannot establish "
                "biological rooting)"
            )
        if not rooted:
            raise TreeInvariantError("species tree is declared unrooted; a rooted tree is required")
    if require_bifurcating_root and len(root_children) != 2:
        raise TreeInvariantError(
            f"species tree root has {len(root_children)} children; a bifurcating root "
            "requires exactly 2 (note: this checks structure only, not biological rooting)"
        )
