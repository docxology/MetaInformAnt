from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

from .amalgkit import (
    AmalgkitParams,
    run_amalgkit,
)


@dataclass(slots=True)
class AmalgkitWorkflowConfig:
    work_dir: Path
    threads: int = 4
    species_list: list[str] = field(default_factory=list)
    # Additional common parameters can be added as needed

    def to_common_params(self) -> AmalgkitParams:
        params: dict[str, Any] = {"threads": self.threads}
        if self.species_list:
            params["species-list"] = list(self.species_list)
        return params


def plan_workflow(config: AmalgkitWorkflowConfig) -> list[tuple[str, AmalgkitParams]]:
    """Return an ordered list of (subcommand, params) representing a full run.

    This does not execute anything; it allows dry inspection and TDD.
    """
    common = config.to_common_params()

    def merge_params(extra: Mapping[str, Any] | None = None) -> AmalgkitParams:
        if not extra:
            return dict(common)
        merged = dict(common)
        merged.update(extra)
        return merged

    steps: list[tuple[str, AmalgkitParams]] = [
        ("metadata", merge_params({})),
        ("integrate", merge_params({})),
        ("config", merge_params({})),
        ("select", merge_params({})),
        ("getfastq", merge_params({})),
        ("quant", merge_params({})),
        ("merge", merge_params({})),
        ("cstmm", merge_params({})),
        ("curate", merge_params({})),
        ("csca", merge_params({})),
        ("sanity", merge_params({})),
    ]
    return steps


def execute_workflow(config: AmalgkitWorkflowConfig, *, check: bool = False) -> list[int]:
    """Execute the full amalgkit workflow in order.

    Returns a list of return codes per step in order.
    """
    steps = plan_workflow(config)
    return_codes: list[int] = []
    for subcommand, params in steps:
        result = run_amalgkit(
            subcommand,
            params,
            work_dir=config.work_dir,
            check=check,
            capture_output=True,
        )
        return_codes.append(result.returncode)
        if check and result.returncode != 0:
            break
    return return_codes


__all__ = [
    "AmalgkitWorkflowConfig",
    "plan_workflow",
    "execute_workflow",
]


