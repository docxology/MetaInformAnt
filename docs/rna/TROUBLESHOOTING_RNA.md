# RNA Pipeline Troubleshooting

Common issues and solutions for running RNA-seq workflows with METAINFORMANT.

**See**: [GETTING_STARTED.md](GETTING_STARTED.md) for setup and usage.

---

## "amalgkit: command not found"

```bash
# Ensure virtual environment is activated
source .venv/bin/activate

# Reinstall if needed with uv
uv pip install --force-reinstall git+https://github.com/kfuku52/amalgkit --python .venv/bin/python3
```

## "fasterq-dump: command not found"

```bash
# Install SRA Toolkit
sudo apt-get install -y sra-toolkit

# Or add to PATH if manually installed
export PATH=$PATH:/path/to/sratoolkit/bin
```

## "Permission denied" errors

```bash
# Ensure output directory is writable
chmod -R u+w output/

# Run script from repository root
cd /path/to/MetaInformAnt
python3 scripts/rna/run_workflow.py ...
```

## External Drive and Filesystem Issues

**See [External Drive Setup Guide](EXTERNAL_DRIVE_SETUP.md) for comprehensive documentation.**

Common issues on external drives (ext6 filesystems):
- Symlink errors when installing packages
- Virtual environment creation failures
- UV cache directory issues

**All issues are automatically handled** by the setup utilities, but see the guide for manual troubleshooting.

## PEP 668 "externally-managed-environment" error

This is expected on modern Linux systems. Always use virtual environments with uv:

```bash
# Create venv with uv
uv venv

# Activate venv
source .venv/bin/activate

# Install packages with uv (no activation needed if using --python flag)
uv pip install -e . --python .venv/bin/python3
```

## Downloads Failing

- Check network connectivity: `ping -c 4 8.8.8.8`
- Verify NCBI_EMAIL is set: `echo $NCBI_EMAIL`
- Try reducing threads: set `threads: 16` in the species config, or `export AK_THREADS=16`
- Check disk space: `df -h /`
- Test ENA connectivity: `wget --spider https://www.ebi.ac.uk/ena/`

## Disk Space Issues

- Clean up partial downloads: `python3 scripts/rna/cleanup_partial_downloads.py --execute`
- Reduce `num_download_workers` in config file for fewer parallel downloads
- Use immediate per-sample processing (default): only one sample's FASTQs exist at a time

## Workflow Not Starting

**Check virtual environment:**
```bash
# Verify venv exists (check both locations)
ls -la .venv/ || ls -la /tmp/metainformant_venv/

# Verify amalgkit installed
source .venv/bin/activate || source /tmp/metainformant_venv/bin/activate
amalgkit --version
```

**Check dependencies:**
```bash
# Verify wget
which wget
wget --version

# Verify kallisto
which kallisto
kallisto version
```

## Processes Stuck or Hung

**Kill and restart:**
```bash
# Find process IDs
ps aux | grep "run_workflow\|amalgkit" | grep -v grep

# Kill specific workflow
kill <PID>

# Or kill all workflows (use with caution)
pkill -f run_workflow

# Restart (will resume from last completed sample)
nohup python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml > logs/species.log 2>&1 &
```

## Quantification Failing

**Check kallisto index:**
```bash
# Verify index exists
ls -lh output/amalgkit/*/work/index/*.idx

# Rebuild if needed (automatic on first run)
# Just restart the workflow
```

**Check FASTQ quality:**
```bash
# List downloaded FASTQs
ls -lh output/amalgkit/*/fastq/*/

# Check if files are complete (not 0 bytes)
find output/amalgkit/*/fastq -name "*.fastq.gz" -size 0
```

## Curate Completes in 0 Seconds

**Cause**: Existing outputs detected, no reprocessing needed

**Solution**: This is normal! To regenerate:
```bash
rm -rf work/curate/*
amalgkit curate --out_dir work --batch_effect_alg no
```

## R Package Errors (Rtsne, amap, vegan, sva, RUVSeq)

**Cause**: System-level gcc compilation issues (documented)

**Solution**: Already handled via patches to `amalgkit/curate.r`:
- Changed `library()` to `require()` for graceful degradation
- Added fallbacks for missing functions
- Core visualizations still work perfectly ✅

**Impact**: Minimal — only advanced optional features affected

---

**See Also**: [GETTING_STARTED.md](GETTING_STARTED.md) | [cloud/TROUBLESHOOTING.md](../cloud/TROUBLESHOOTING.md) | [Steps Index](amalgkit/steps/README.md)
