# Running End-to-End Amalgkit Workflows for All Species

**Quick Answer**: Use `run_multi_species.py` for all species with configurable threads.

## Prerequisites

### Required Setup

Before running workflows, ensure you have:

1. **Virtual Environment**: Created and activated using `uv` (recommended)
   ```bash
   # Install uv if not already installed
   curl -LsSf https://astral.sh/uv/install.sh | sh
   
   # Create virtual environment (uses /tmp/metainformant_venv on ext6 filesystems)
   uv venv .venv
   
   # Activate it
   source .venv/bin/activate
   # Or if using alternative location: source /tmp/metainformant_venv/bin/activate
   
   # Install dependencies
   uv pip install -e .
   uv pip install git+https://github.com/kfuku52/amalgkit
   ```
   
   **Note**: On ext6 filesystems (which don't support symlinks), the setup automatically uses `/tmp/metainformant_venv` if `.venv` creation fails. Scripts automatically discover and use the venv location.

2. **Verify Installation**:
   ```bash
   # Check amalgkit is installed
   amalgkit --version
   
   # Check other tools
   which wget kallisto fastp seqkit
   ```

3. **Environment Variables** (optional but recommended):
   ```bash
   # Set NCBI email for API calls
   export NCBI_EMAIL="your@email.com"
   
   # Set threads globally (optional)
   export AK_THREADS=12
   ```

**Note**: Scripts automatically discover and activate virtual environments (`.venv` or `/tmp/metainformant_venv`). The venv must be set up first using `uv venv` and `uv pip install`.

## Quick Start

### Method 1: All Species with Configurable Threads (Recommended)

**For end-to-end workflows (metadata → sanity) on all species:**

```bash
# Ensure virtual environment is set up (see Prerequisites above)

# Run all species with threads from config files (default: 10 per species)
python3 scripts/rna/run_multi_species.py

# To customize threads, use environment variable:
export AK_THREADS=12
python3 scripts/rna/run_multi_species.py
```

**What it does:**
- ✅ Auto-discovers all `config/amalgkit/amalgkit_*.yaml` files
- ✅ Runs full end-to-end workflow (metadata → select → getfastq → quant → merge → curate → sanity)
- ✅ Processes species sequentially (one at a time)
- ✅ Uses threads from config files (default: 10) or `AK_THREADS` environment variable
- ✅ Auto-activates virtual environment
- ✅ Cross-species analysis (CSTMM, CSCA) after all species complete

### Method 2: Loop Over All Configs (ENA Workflow)

**For maximum reliability with ENA downloads:**

```bash
# Run all species in parallel (background processes)
for config in config/amalgkit/amalgkit_*.yaml; do
  # Skip template
  [[ "$config" == *template* ]] && continue
  
  species=$(basename "$config" .yaml | sed 's/amalgkit_//')
  nohup python3 scripts/rna/workflow_ena_integrated.py \
    --config "$config" \
    --batch-size 12 \
    --threads 12 \
    > output/workflow_${species}_$(date +%Y%m%d_%H%M%S).log 2>&1 &
done

# Or run sequentially
for config in config/amalgkit/amalgkit_*.yaml; do
  [[ "$config" == *template* ]] && continue
  python3 scripts/rna/workflow_ena_integrated.py \
    --config "$config" \
    --batch-size 12 \
    --threads 12
done
```

**What it does:**
- ✅ Uses ENA direct downloads (100% reliability)
- ✅ Configurable threads via `--threads` argument
- ✅ Batched processing (12 samples at a time)
- ✅ Auto-cleanup (FASTQs deleted after quantification)
- ⚠️ Requires manual loop for multiple species

## Thread Configuration

### Option 1: Environment Variable (Recommended)

```bash
# Set threads globally (affects all workflows)
export AK_THREADS=12
python3 scripts/rna/run_multi_species.py
```

### Option 2: Edit Config Files

Edit `threads:` in each config file:
```yaml
# config/amalgkit/amalgkit_cfloridanus.yaml
threads: 12  # Change from default 10 to 12
```

### Option 3: Command-Line Override (ENA Workflow Only)

```bash
# workflow_ena_integrated.py supports --threads argument
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --threads 12
```

## Complete Examples

### Example 1: All Species with 12 Threads (SRA Workflow)

```bash
# 1. Ensure virtual environment is set up
if [ ! -d ".venv" ] && [ ! -d "/tmp/metainformant_venv" ]; then
    echo "Creating virtual environment with uv..."
    uv venv .venv || uv venv /tmp/metainformant_venv
    source .venv/bin/activate || source /tmp/metainformant_venv/bin/activate
    uv pip install -e .
    uv pip install git+https://github.com/kfuku52/amalgkit
else
    # Activate existing venv (scripts auto-discover location)
    [ -d ".venv" ] && source .venv/bin/activate || source /tmp/metainformant_venv/bin/activate
fi

# 2. Set threads via environment variable
export AK_THREADS=12

# 3. Run all species
python3 scripts/rna/run_multi_species.py
```

### Example 2: All Species with 12 Threads (ENA Workflow - Parallel)

```bash
# Create a helper script
cat > run_all_species_ena.sh << 'EOF'
#!/bin/bash
for config in config/amalgkit/amalgkit_*.yaml; do
  [[ "$config" == *template* ]] && continue
  [[ "$config" == *test* ]] && continue
  
  species=$(basename "$config" .yaml | sed 's/amalgkit_//')
  echo "Starting $species..."
  
  nohup python3 scripts/rna/workflow_ena_integrated.py \
    --config "$config" \
    --batch-size 12 \
    --threads 12 \
    > output/workflow_${species}_$(date +%Y%m%d_%H%M%S).log 2>&1 &
done

echo "All workflows started. Check logs in output/"
EOF

chmod +x run_all_species_ena.sh
./run_all_species_ena.sh
```

### Example 3: All Species with 12 Threads (ENA Workflow - Sequential)

```bash
# Process one species at a time
for config in config/amalgkit/amalgkit_*.yaml; do
  [[ "$config" == *template* ]] && continue
  [[ "$config" == *test* ]] && continue
  
  species=$(basename "$config" .yaml | sed 's/amalgkit_//')
  echo "Processing $species..."
  
  python3 scripts/rna/workflow_ena_integrated.py \
    --config "$config" \
    --batch-size 12 \
    --threads 12
done
```

## Comparison

| Method | Script | Threads Config | End-to-End | All Species | Reliability |
|--------|--------|----------------|------------|-------------|-------------|
| **Method 1** | `run_multi_species.py` | Config file or `AK_THREADS` | ✅ Yes | ✅ Auto | ~0% (SRA) |
| **Method 2** | `workflow_ena_integrated.py` (loop) | `--threads` argument | ✅ Yes | ⚠️ Manual loop | 100% (ENA) |

## Recommended Approach

**For production (all species, configurable threads):**

1. **If you need reliability**: Use Method 2 (ENA workflow) with a loop
2. **If you need simplicity**: Use Method 1 (`run_multi_species.py`) with `AK_THREADS=12`

**Example:**
```bash
# 1. Setup (one-time, if not already done)
uv venv .venv  # or /tmp/metainformant_venv on ext6 filesystems
source .venv/bin/activate  # or /tmp/metainformant_venv/bin/activate
uv pip install -e .
uv pip install git+https://github.com/kfuku52/amalgkit

# 2. Run all species with configurable threads
# Note: Scripts automatically discover venv location (.venv or /tmp/metainformant_venv)
export AK_THREADS=12
python3 scripts/rna/run_all_species_parallel.py  # Recommended: parallel execution
# Or: python3 scripts/rna/run_multi_species.py  # Sequential execution
```

## Monitoring

```bash
# Check status (unified orchestrator)
python3 scripts/rna/orchestrate_workflows.py --status

# Real-time monitoring
python3 scripts/rna/orchestrate_workflows.py --monitor

# Detailed status with categories
python3 scripts/rna/orchestrate_workflows.py --status --detailed

# Check running processes
ps aux | grep -E "(workflow_ena|run_all_species_parallel|run_multi_species)" | grep -v grep

# View logs
tail -f output/parallel_all_species_*.log
tail -f output/workflow_*.log
```

## See Also

- **[ORCHESTRATION/README.md](ORCHESTRATION/README.md)**: Overview of all orchestrators
- **[MULTI_SPECIES_QUICK_START.md](MULTI_SPECIES_QUICK_START.md)**: Detailed production guide
- **[CONFIGURATION.md](CONFIGURATION.md)**: Configuration management

