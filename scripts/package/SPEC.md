# SPEC: Package Scripts

Utility scripts for MetaInformAnt package management, dependency handles, and repository-wide test orchestration.

## Core Workflows

- `run_all_scripts.py`: The master orchestrator for discovering and executing all module-level thin wrappers.
- `setup_environment.sh`: Bootstrapping script for setting up `uv` and biological dependencies.
- `check_module_health.py`: Validates documentation compliance and basic import health for all modules.

## Standards

- **Automation**: Scripts in this directory should be idempotent where possible.
- **Reporting**: Comprehensive logging of successes and failures is mandatory for the master orchestrator.
