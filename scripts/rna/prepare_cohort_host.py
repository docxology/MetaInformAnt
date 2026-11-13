#!/usr/bin/env python3
"""Prepare a fresh data root for frozen-cohort quantification on a new host.

Thin entry point; see ``metainformant.rna.engine.cohort_host``. Places the
cohort's exact kallisto index bytes where the species configuration resolves
them and writes the reference manifests cohort discovery and imports require.
The default is a read-only plan.

Example::

    python3 scripts/rna/prepare_cohort_host.py --data-root ROOT \\
        --config-dir projects/hymenoptera_amalgkit/config/amalgkit \\
        --cohort ROOT/cohort/hymenoptera27-20261005 --species monomorium_pharaonis --execute
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.rna.engine.cohort_host import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
