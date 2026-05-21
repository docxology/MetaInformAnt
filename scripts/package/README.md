# PACKAGE

## Overview
Command-line helpers for Package workflows. Scripts should remain thin wrappers around `src/metainformant/` implementations and be run from the repository root with `uv`.

## Contents
- [_common.sh](_common.sh)
- [build.sh](build.sh)
- [build_utils.sh](build_utils.sh)
- [fix_tmp_space.sh](fix_tmp_space.sh)
- [generate_cursor_skills.py](generate_cursor_skills.py)
- [generate_custom_summary.py](generate_custom_summary.py)
- [install_linux_deps.sh](install_linux_deps.sh)
- [patch_amalgkit.py](patch_amalgkit.py)
- [release.sh](release.sh)
- [run_tests.sh](run_tests.sh)
- [setup.sh](setup.sh)
- [setup_uv.sh](setup_uv.sh)
- [summarize_status.py](summarize_status.py)
- [test.sh](test.sh)
- [uv_dev_setup.sh](uv_dev_setup.sh)
- [uv_docs.sh](uv_docs.sh)
- [uv_profile.sh](uv_profile.sh)
- [uv_quality.sh](uv_quality.sh)
- [uv_test.sh](uv_test.sh)
- [uv_test_optimized.sh](uv_test_optimized.sh)
- [uv_test_setup.sh](uv_test_setup.sh)
- [validate_build.sh](validate_build.sh)
- [verify.sh](verify.sh)
- [verify_test_deps.sh](verify_test_deps.sh)
- [verify_uv_setup.sh](verify_uv_setup.sh)

## Usage
```bash
uv run python scripts/package/generate_cursor_skills.py --help
```
