"""Public Amalgkit acquisition API; scientific and scheduling methods live in the engine."""
from __future__ import annotations
from metainformant.rna.engine.acquisition_allocation import AcquisitionAllocation, allocate_tasks, write_allocation
from metainformant.rna.engine.acquisition_estimates import ThroughputEvidence, LaneCosts, LaneEstimate, CampaignEstimate, estimate_lane, combine_estimates
from metainformant.rna.engine.acquisition_snapshot import create_campaign_manifest, stage_worker_inputs
from metainformant.rna.engine.acquisition_worker import run_manifest
from metainformant.rna.engine.acquisition_cli import main

__all__ = ["AcquisitionAllocation", "allocate_tasks", "write_allocation", "ThroughputEvidence", "LaneCosts",
           "LaneEstimate", "CampaignEstimate", "estimate_lane", "combine_estimates", "create_campaign_manifest", "stage_worker_inputs", "run_manifest"]

if __name__ == "__main__":
    raise SystemExit(main())
