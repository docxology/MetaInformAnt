from __future__ import annotations

import argparse
import json
from pathlib import Path

from metainformant.rna.engine.campaign_observation import capture_campaign_observation
from metainformant.rna.engine.campaign_observation_report import (
    write_observation_bundle,
)


def main(
    argv: list[str] | None = None, *, default_title: str = "RNA-seq acquisition"
) -> int:
    parser = argparse.ArgumentParser(
        description="Capture a read-only operational observation and isolated working supplement."
    )
    parser.add_argument("--operations-root", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--transfer-observation", type=Path)
    parser.add_argument("--task-selection", type=Path, required=True)
    parser.add_argument("--journal", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-rate-age-seconds", type=float, default=300)
    parser.add_argument("--title", default=default_title)
    parser.add_argument(
        "--no-supplement-pdf",
        action="store_true",
        help="Keep Markdown and figure exports; skip pandoc supplement.",
    )
    args = parser.parse_args(argv)
    observation = capture_campaign_observation(
        args.operations_root,
        args.data_root,
        selection_path=args.task_selection,
        journal_path=args.journal,
        manifest_path=args.manifest,
        transfer_observation_path=args.transfer_observation,
        max_rate_age_seconds=args.max_rate_age_seconds,
    )
    output = write_observation_bundle(
        observation,
        args.output_dir,
        render_pdf=not args.no_supplement_pdf,
        title=args.title,
    )
    print(
        json.dumps(
            {
                "output_dir": str(output),
                "counts": observation.payload["counts"],
                "scientific_release_ready": False,
            },
            indent=2,
        )
    )
    return 0
