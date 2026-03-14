# Bioinformatics Project Template

This repository serves as a **standard template** for organizing bioinformatics projects.
It ensures reproducibility, clear separation of data and code, and organized tracking of outputs and execution parameters.

## Quickstart

This boilerplate provides everything you need to start a data-processing pipeline.

1. **Environment Setup**
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
    ```

2. **Configuration**
   Open `config/default.yaml` to define your input bounds, thresholds, and logging formatting without hardcoding variables into scripts.

3. **Execution**
   A standard entrypoint script is provided:
    ```bash
    python scripts/01_process_data.py --config config/default.yaml
    ```
   This script will read from `data/raw/` and establish its log tracking in `logs/`.

## Directory Structure

Strict adherence to the directory structure preserves project integrity:

- `data/`
  - `raw/`: **Immutable raw data.** Treat everything in this folder as absolutely read-only.
  - `processed/`: Intermediate cleansed and processed data.
- `scripts/`: Production source code for data processing and downstream analysis.
- `notebooks/`: Jupyter notebooks (`.ipynb`) restricted purely to exploratory data analysis and visualization.
- `results/`: Formatted outputs (e.g., final TSV tables, figures) from scripts and notebooks.
- `logs/`: High-resolution execution logs generated directly by scripts to track state over time.
- `config/`: Pipeline configuration files (YAML, JSON) specifying hyper-parameters and path routing.

## Principles of Reproducibility

1. **No Hardcoding**: All paths and thresholds should route dynamically through `config/` files.
2. **Immutability**: Anything inside `data/raw` must never be modified by a script. Generate altered states into `data/processed`.
3. **Traceability**: All scripts should utilize the standard Python `logging` library mapping `.log` files securely into `logs/`.
