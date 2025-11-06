#!/bin/bash
# Complete workflow sequence for Pogonomyrmex barbatus
# Run these commands sequentially

set -e  # Exit on error

CONFIG="config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml"

echo "============================================================================"
echo "POGONOMYRMEX BARBATUS WORKFLOW - SEQUENTIAL EXECUTION"
echo "============================================================================"
echo ""

# Step 1: Check Environment
echo "Step 1: Checking environment..."
python3 scripts/rna/check_environment.py
echo ""

# Step 2: Verify Genome Setup
echo "Step 2: Verifying genome setup..."
python3 scripts/rna/setup_genome.py --config "$CONFIG" --verify-only
echo ""

# Step 3: Check Current Status
echo "Step 3: Checking current workflow status..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --status --detailed
echo ""

# Step 4: Metadata (discover samples from SRA)
echo "Step 4: Running metadata step (discover samples from SRA)..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --steps metadata
echo ""

# Step 5: Integrate (prepare metadata - will skip if no FASTQ files yet)
echo "Step 5: Running integrate step..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --steps integrate
echo ""

# Step 6: Config (generate config files)
echo "Step 6: Running config step..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --steps config
echo ""

# Step 7: Select (filter samples)
echo "Step 7: Running select step (filter samples)..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --steps select
echo ""

# Step 8: Getfastq (download FASTQ files - THIS IS THE LONGEST STEP)
echo "Step 8: Running getfastq step (download FASTQ files)..."
echo "WARNING: This step can take hours or days depending on number of samples"
python3 scripts/rna/run_workflow.py --config "$CONFIG" --steps getfastq
echo ""

# Step 9: Integrate (now that FASTQ files exist)
echo "Step 9: Running integrate step again (now with FASTQ files)..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --steps integrate
echo ""

# Step 10: Quant (quantify with kallisto)
echo "Step 10: Running quant step (quantify with kallisto)..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --steps quant
echo ""

# Step 11: Merge (merge quantification results)
echo "Step 11: Running merge step (merge quantification results)..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --steps merge
echo ""

# Step 12: Curate (quality control)
echo "Step 12: Running curate step (quality control)..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --steps curate
echo ""

# Step 13: Sanity (validate outputs)
echo "Step 13: Running sanity step (validate outputs)..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --steps sanity
echo ""

# Step 14: Final Status Check
echo "Step 14: Final status check..."
python3 scripts/rna/run_workflow.py --config "$CONFIG" --status --detailed
echo ""

echo "============================================================================"
echo "WORKFLOW COMPLETE!"
echo "============================================================================"

