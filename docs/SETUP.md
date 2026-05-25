# Installation & Setup Guide

Complete environment setup for METAINFORMANT across platforms and use cases.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Platform-Specific Setup](#platform-specific-setup)
  - [Linux (Ubuntu/Debian)](#linux-ubuntudebian)
  - [macOS](#macos)
  - [Windows + WSL2](#windows--wsl2)
  - [External Drive (exFAT/FAT32)](#external-drive-exfatfat32)
- [Installation Methods](#installation-methods)
  - [Automated Setup (Recommended)](#automated-setup-recommended)
  - [Manual Setup](#manual-setup)
  - [Development Mode](#development-mode)
- [Post-Installation](#post-installation)
  - [Verify Installation](#verify-installation)
  - [Download Test Data](#download-test-data)
  - [Run Demo](#run-demo)
- [Troubleshooting](#troubleshooting)
- [Next Steps](#next-steps)

---

## Prerequisites

| Tool | Purpose | Install Command |
|------|---------|-----------------|
| **Python 3.11+** | Runtime | `sudo apt install python3.11` (Ubuntu) |
| **`uv`** | Package manager (required) | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| **Git** | Version control | `sudo apt install git` |
| **gcloud CLI** *(optional)* | GCP cloud deployment | `bash scripts/cloud/install_gcloud.sh` |

### Why `uv`?

- Lightning-fast dependency resolution (Rust-based)
- Composable virtual environments
- Lockfile generation for reproducible installs
- Works on FAT filesystems via cache redirection
- Drop-in replacement for `pip` + `venv`

**Never use `pip install` directly** — see [REAL_IMPLEMENTATION_POLICY.md](REAL_IMPLEMENTATION_POLICY.md).

---

## Platform-Specific Setup

### Linux (Ubuntu/Debian)

```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.bashrc  # or restart shell

# Clone and setup
git clone https://github.com/docxology/metainformant.git
cd metainformant
bash scripts/package/setup.sh

# Activate
source .venv/bin/activate

# Verify
python -c "import metainformant; print(metainformant.__version__)"
```

### macOS

```bash
# Install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install uv and Python
brew install uv python@3.11

# Clone and setup
git clone https://github.com/docxology/metainformant.git
cd metainformant
bash scripts/package/setup.sh

# Activate
source .venv/bin/activate
```

### Windows + WSL2

WSL2 required — native Windows Python fails due to symlink limitations.

```powershell
# PowerShell (admin): Install WSL Ubuntu
wsl --install -d Ubuntu

# In WSL Ubuntu terminal:
sudo apt update
curl -LsSf https://astral.sh/uv/install.sh | sh
git clone https://github.com/docxology/metainformant.git
cd metainformant
bash scripts/package/setup.sh
source .venv/bin/activate
```

**Note:** Never install in `/mnt/c/...` (Windows filesystem). Use WSL `~/` path.

---

## Installation Methods

### Automated Setup (Recommended)

`bash scripts/package/setup.sh` handles:
- Python/uv detection
- venv creation
- Dependency installation
- Test data download (optional)
- Smoke test suite

### Manual Setup

```bash
uv venv
source .venv/bin/activate
uv pip install -e ".[dev,all]"
```

### Development Mode

```bash
uv pip install -e ".[dev,test,docs]"
pre-commit install  # Enable git hooks
```

---

## Post-Installation

### Verify Installation

```bash
# Check Python
python --version  # 3.11+

# Check package
python -c "import metainformant; print(metainformant.__version__)"

# Run test suite (optional but recommended)
pytest tests/ -v --tb=short
```

### Download Test Data

```bash
bash scripts/data/download_test_datasets.sh
```

### Run Demo

```bash
python3 scripts/core/run_demo.py
```

---

## Troubleshooting

| Issue | Cause | Fix |
|-------|-------|-----|
| `uv: command not found` | Not in PATH | `source ~/.bashrc` or reinstall |
| `Python 3.9 detected` | System Python takes precedence | `export PATH="/opt/homebrew/bin:$PATH"` (macOS) |
| venv creation fails on FAT drive | Symlink not supported | Use `/tmp` venv: `export UV_CACHE_DIR=/tmp/uv-cache` |
| ImportError: No module named X | Dependencies not installed | `uv pip install -e ".[all]"` |
| `Permission denied` | Wrong Python interpreter | Check `which python` is from venv |

**Logs:** Check `~/.hermes/logs/agent.log` if using Hermes agent.

---

## Next Steps

1. Read [TUTORIALS.md](TUTORIALS.md) for hands-on examples
2. Pick your [module](../README.md#choosing-the-right-module)
3. Join [community](../CONTRIBUTING.md#community)

---

**Related:** [UV_SETUP.md](UV_SETUP.md) | [EXTERNAL_DRIVE_SETUP](rna/EXTERNAL_DRIVE_SETUP.md) | [real-implementation policy](REAL_IMPLEMENTATION_POLICY.md)
