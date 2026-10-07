"""Public Amalgkit acquisition API; scientific and scheduling methods live in the engine."""

from __future__ import annotations

from metainformant.rna.engine.acquisition_allocation import AcquisitionAllocation, allocate_tasks, write_allocation
from metainformant.rna.engine.acquisition_cli import main
from metainformant.rna.engine.acquisition_estimates import (
    CampaignEstimate,
    LaneCosts,
    LaneEstimate,
    ThroughputEvidence,
    combine_estimates,
    estimate_lane,
)
from metainformant.rna.engine.acquisition_snapshot import create_campaign_manifest, stage_worker_inputs
from metainformant.rna.engine.acquisition_worker import run_manifest

__all__ = [
    "AcquisitionAllocation",
    "allocate_tasks",
    "write_allocation",
    "ThroughputEvidence",
    "LaneCosts",
    "LaneEstimate",
    "CampaignEstimate",
    "estimate_lane",
    "combine_estimates",
    "create_campaign_manifest",
    "stage_worker_inputs",
    "run_manifest",
]

if __name__ == "__main__":
    raise SystemExit(main())
