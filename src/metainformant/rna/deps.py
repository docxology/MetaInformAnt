from __future__ import annotations

"""Lightweight external dependency checks for RNA workflow steps.

This module reports whether non-Python, external CLIs that particular
`amalgkit` steps rely upon are present on PATH. The checks are intentionally
conservative, fast, and offline.

Policy:
- Return ok=True when at least one acceptable alternative is present.
- Do not perform network calls.
- Keep step-specific knowledge localized here to avoid scattering checks.
"""

from dataclasses import dataclass
from shutil import which
from typing import Iterable


@dataclass(frozen=True)
class StepDependencyStatus:
    """Status of step dependency checking.
    
    Attributes:
        ok: True if all required dependencies are present
        missing: List of missing dependency names
        present: List of present dependency names
        note: Optional note about status
    """
    ok: bool
    missing: list[str]
    present: list[str]
    note: str = ""


def _any_present(candidates: Iterable[str]) -> tuple[bool, list[str], list[str]]:
    present: list[str] = []
    missing: list[str] = []
    any_ok = False
    for name in candidates:
        if which(name) is not None:
            present.append(name)
            any_ok = True
        else:
            missing.append(name)
    return any_ok, missing, present


def check_step_dependencies(step: str) -> StepDependencyStatus:
    """Return external CLI availability for a given step.

    Steps and their checks:
    - metadata, integrate, config, select, merge, sanity: rely on `amalgkit` which
      is checked separately in the workflow preflight; no extra hard deps here.
    - getfastq: require at least one of parallel-fastq-dump or sra-tools utilities
      (prefetch/fastq-dump/fasterq-dump).
    - quant: require at least one of salmon or kallisto.
    - cstmm/curate/csca: require R (Rscript available) for plotting/stats.
    """
    step = step.strip().lower()

    if step in {"metadata", "config", "select", "merge", "sanity"}:
        return StepDependencyStatus(ok=True, missing=[], present=[], note="no extra deps")

    if step == "integrate":
        ok, missing, present = _any_present(["seqkit"])
        return StepDependencyStatus(
            ok=ok,
            missing=[] if ok else missing,
            present=present,
            note=("need seqkit for integrate step"),
        )

    if step == "getfastq":
        # Require either sratoolkit tools present, or both parallel-fastq-dump and an sra-dump tool
        have_pfd, _, present_pfd = _any_present(["parallel-fastq-dump"])
        have_sra_tool, _, present_sra = _any_present(["prefetch", "fastq-dump", "fasterq-dump"])  # any one suffices
        ok = have_sra_tool or (have_pfd and have_sra_tool)
        missing: list[str] = [] if ok else ["sratoolkit (prefetch/fastq-dump/fasterq-dump)"]
        present = present_pfd + present_sra
        return StepDependencyStatus(
            ok=ok,
            missing=missing,
            present=present,
            note=("need sratoolkit executables; parallel-fastq-dump alone is insufficient"),
        )

    if step == "quant":
        ok, missing, present = _any_present(["salmon", "kallisto"])
        return StepDependencyStatus(
            ok=ok,
            missing=[] if ok else missing,
            present=present,
            note=("need salmon or kallisto"),
        )

    if step in {"cstmm", "curate", "csca"}:
        ok, missing, present = _any_present(["Rscript", "R"])
        return StepDependencyStatus(
            ok=ok,
            missing=[] if ok else missing,
            present=present,
            note=("need R (Rscript) for downstream analyses and plotting"),
        )

    # Unknown step: be permissive
    return StepDependencyStatus(ok=True, missing=[], present=[], note="unknown step")


__all__ = ["StepDependencyStatus", "check_step_dependencies"]
