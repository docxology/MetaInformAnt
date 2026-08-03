#!/usr/bin/env python3
"""Verify the installed Amalgkit runtime without modifying site-packages."""

from __future__ import annotations

import inspect
import json
import sys
from typing import Any

from metainformant.rna.amalgkit import (
    AMALGKIT_INSTALL_SPEC,
    REQUIRED_AMALGKIT_VERSION,
    validate_amalgkit_version,
)


def inspect_runtime() -> dict[str, Any]:
    """Return read-only capability facts for the pinned Amalgkit release."""

    import amalgkit.getfastq as getfastq_module
    import amalgkit.subprocess_utils as subprocess_utils

    execute = getattr(getfastq_module, "execute_fasterq_dump_command", None)
    runner = getattr(subprocess_utils, "run_logged_command", None)
    execute_source = inspect.getsource(execute) if execute is not None else ""
    runner_signature = str(inspect.signature(runner)) if runner is not None else None

    return {
        "getfastq_module": getattr(getfastq_module, "__file__", None),
        "subprocess_utils_module": getattr(subprocess_utils, "__file__", None),
        "execute_fasterq_dump_command": execute is not None,
        "getfastq_uses_logged_runner": "run_logged_command" in execute_source,
        "logged_runner_signature": runner_signature,
        "logged_runner_accepts_runner": runner is not None and "runner" in inspect.signature(runner).parameters,
    }


def main() -> int:
    version_ok, version_message = validate_amalgkit_version(REQUIRED_AMALGKIT_VERSION, exact=True)
    payload: dict[str, Any] = {
        "required_version": REQUIRED_AMALGKIT_VERSION,
        "install_spec": AMALGKIT_INSTALL_SPEC,
        "version_ok": version_ok,
        "version_message": version_message,
    }
    if version_ok:
        try:
            payload["capabilities"] = inspect_runtime()
        except (ImportError, OSError, TypeError, ValueError) as exc:
            payload["capability_error"] = str(exc)

    print(json.dumps(payload, indent=2, sort_keys=True))
    capabilities = payload.get("capabilities", {})
    capability_ok = all(
        capabilities.get(key) is True
        for key in (
            "execute_fasterq_dump_command",
            "getfastq_uses_logged_runner",
            "logged_runner_accepts_runner",
        )
    )
    return 0 if version_ok and capability_ok else 1


if __name__ == "__main__":
    sys.exit(main())
