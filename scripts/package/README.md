# PACKAGE

## Overview
Command-line helpers for Package workflows. Scripts should remain thin wrappers around `src/metainformant/` implementations and be run from the repository root with `uv`.

## Contents
- [_common.sh](_common.sh)
- [build.sh](build.sh)
- [build_utils.sh](build_utils.sh)
- [fix_tmp_space.sh](fix_tmp_space.sh)
- [generate_custom_summary.py](generate_custom_summary.py)
- [install_linux_deps.sh](install_linux_deps.sh)
- [verify_amalgkit_runtime.py](verify_amalgkit_runtime.py)
- [release.sh](release.sh)
- [setup.sh](setup.sh)
- [summarize_status.py](summarize_status.py)
- [test.sh](test.sh)
- [uv_dev_setup.sh](uv_dev_setup.sh)
- [uv_docs.sh](uv_docs.sh)
- [uv_profile.sh](uv_profile.sh)
- [uv_quality.sh](uv_quality.sh)
- [validate_build.sh](validate_build.sh)
- [verify.sh](verify.sh)
- [verify_test_deps.sh](verify_test_deps.sh)
- [verify_uv_setup.sh](verify_uv_setup.sh)

## Usage
```bash
bash scripts/package/test.sh                 # fast test lane
bash scripts/package/test.sh --mode coverage # full suite with the coverage floor
bash scripts/package/uv_quality.sh           # format, lint and type checks
```
