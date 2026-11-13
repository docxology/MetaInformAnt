#!/usr/bin/env python3
"""Export verified, reusable quantification outputs to a portable snapshot.

Thin entry point; see ``metainformant.rna.engine.quant_snapshot`` for the
evidence rules, the allowlist and the idempotency contract. The default is a
read-only plan.

Example::

    python3 scripts/rna/export_quant_snapshot.py --data-root ROOT \\
        --config-dir projects/hymenoptera_amalgkit/config/amalgkit \\
        --cohort ROOT/cohort/hymenoptera27-20261005 \\
        --dest projects/hymenoptera_amalgkit/quant/completed/hymenoptera27-20261005 \\
        --git-tree projects/hymenoptera_amalgkit:HEAD \\
        --inventory-dir projects/hymenoptera_amalgkit/quant/inventory/hymenoptera27-20261005
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.rna.engine.quant_snapshot import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
