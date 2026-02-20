---
description: Replicate the MetaInformAnt environment on a new machine
---
# Move Environment to a New Computer

This workflow is specific to replicating the `metainformant` setup onto a new machine, particularly when using a portable drive (like `blue/`).

Because we use `uv`, the process is entirely deterministic, although special consideration must be given to FAT32 or exFAT filesystems which do not support symlinks.

// turbo-all

1. **Clone the Repository** (If not already present on the drive)

   ```bash
   git clone https://github.com/docxology/metainformant.git
   cd metainformant
   ```

2. **Run Initialization Script**

   ```bash
   bash scripts/package/setup.sh
   ```

   *Note: This script automatically detects if the drive uses a FAT filesystem (exFAT/FAT32) and routes the virtual environment to `/tmp/metainformant_venv` and cache to `/tmp/uv-cache` to seamlessly bypass symlink limitations.*

3. **Activate the Environment**
   If on standard Ext4/APFS/NTFS:

   ```bash
   source .venv/bin/activate
   ```

   If on FAT/exFAT:

   ```bash
   source /tmp/metainformant_venv/bin/activate
   ```

4. **Verify Installation**

   ```bash
   python -c "import metainformant; print('Installation Success:', metainformant.__version__)"
   uv run pytest tests/ -v
   ```
