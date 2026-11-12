#!/usr/bin/env python3
"""Reclaim raw inputs only for verified, reusable quantification outputs.

Thin entry point; see ``metainformant.rna.engine.raw_reclamation`` for the
scopes, the evidence rules and the locking contract.  The default is a
read-only dry run.

Example::

    python3 scripts/rna/reclaim_quantified_raw.py --data-root ROOT \\
        --cohort ROOT/cohort/hymenoptera27-20261005 --include-staging --include-sra-cache \\
        --include-excluded --manifest /tmp/reclaim-plan.json
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.rna.engine.raw_reclamation import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
