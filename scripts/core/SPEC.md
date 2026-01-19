# SPEC: Core Scripts

Infrastructure utility scripts for system maintenance and diagnostic reporting.

## Key Scripts

- `check_system_deps.py`: Verifies the presence of all biological and system binaries (R, GATK, etc.).
- `generate_diagnostic_report.py`: Aggregates logs and results for troubleshooting.

## Standards

- **Portability**: Must remain platform-agnostic (macOS/Linux).
- **Automation**: Designed for execution within CI/CD or deployment pipelines.
