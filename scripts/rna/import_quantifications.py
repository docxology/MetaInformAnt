#!/usr/bin/env python3
"""Import verified quantifications from other data roots into a campaign data root.

Dry run (default) writes a JSON-lines decision manifest to ``--manifest`` and
changes nothing; ``--execute`` copies verified samples under the data-root
maintenance lock and never replaces an existing sample directory. See
``metainformant.rna.engine.quant_import`` for the decision rules.

Usage:
    python3 scripts/rna/import_quantifications.py --data-root ROOT --cohort COHORT_DIR \\
        --config-dir projects/hymenoptera_amalgkit/config/amalgkit \\
        --source october=/path/to/october_root --source cloud=/path/to/completed_quant:flat \\
        --receipt-store oct=/path/to/october_root/locked-quant:local-20261019 \\
        --manifest /tmp/import-plan.jsonl
    ... --execute [--adopt-db ROOT/pipeline_progress.db]
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.rna.engine.quant_import import main

if __name__ == "__main__":
    raise SystemExit(main())
