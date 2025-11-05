#!/usr/bin/env bash
# Example: Discover and run RNA-seq workflows for all ant species
#
# This script demonstrates the complete workflow from discovery to execution

set -e  # Exit on error

echo "================================================================================"
echo "ANT SPECIES DISCOVERY AND WORKFLOW EXECUTION EXAMPLE"
echo "================================================================================"
echo ""

# Step 1: Prerequisites check
echo "Step 1: Checking prerequisites..."
if ! python3 -c "import Bio; from ncbi.datasets import GenomeApi" 2>/dev/null; then
    echo "❌ Missing dependencies. Installing..."
    pip install biopython ncbi-datasets-pylib
fi

if [ -z "$NCBI_EMAIL" ]; then
    echo "⚠️  NCBI_EMAIL not set. Using placeholder..."
    export NCBI_EMAIL="your.email@example.com"
    echo "   Please set: export NCBI_EMAIL=\"your@email.com\""
fi

echo "✅ Prerequisites ready"
echo ""

# Step 2: Discover ant species with RNA-seq data
echo "Step 2: Discovering ant species with RNA-seq data..."
echo "   This may take 5-15 minutes..."
python3 scripts/rna/discover_ant_species_with_rnaseq.py \
    --output-dir output/ant_discovery \
    --min-samples 10

echo ""
echo "✅ Discovery complete!"
echo ""

# Step 3: Review discovery report
echo "Step 3: Discovery report summary..."
echo ""
cat output/ant_discovery/DISCOVERY_REPORT.md | head -50
echo ""
echo "   Full report: output/ant_discovery/DISCOVERY_REPORT.md"
echo ""

# Step 4: Deploy high-priority configurations
echo "Step 4: Deploying configurations for high-priority species..."

# Create backup of existing configs
if [ -d "config/amalgkit" ]; then
    echo "   Backing up existing configs..."
    mkdir -p config/amalgkit_backup_$(date +%Y%m%d)
    cp config/amalgkit/*.yaml config/amalgkit_backup_$(date +%Y%m%d)/ 2>/dev/null || true
fi

# Deploy configs for species with chromosome-level genomes and many samples
echo "   Deploying chromosome-level genomes with ≥10 samples..."

# Parse the discovery JSON to find high-priority species
python3 << 'PYEOF'
import json
from pathlib import Path
import shutil

# Load discovery data
with open('output/ant_discovery/ant_species_rnaseq_data.json') as f:
    species_data = json.load(f)

# Find high-priority species
high_priority = []
for species_name, data in species_data.items():
    if data.get('genome'):
        genome = data['genome']
        # Chromosome-level with ≥10 samples
        if genome['level'] == 'Chromosome' and data['sample_count'] >= 10:
            high_priority.append({
                'name': species_name,
                'slug': species_name.replace(' ', '_').lower(),
                'samples': data['sample_count'],
                'accession': genome['accession']
            })

# Sort by sample count
high_priority.sort(key=lambda x: x['samples'], reverse=True)

# Deploy configs
config_dir = Path('config/amalgkit')
config_dir.mkdir(parents=True, exist_ok=True)

print(f"\nDeploying {len(high_priority)} high-priority species:")
for species in high_priority:
    src = Path(f"output/ant_discovery/configs/amalgkit_{species['slug']}.yaml")
    dst = config_dir / f"amalgkit_{species['slug']}.yaml"
    
    if src.exists():
        shutil.copy2(src, dst)
        print(f"  ✅ {species['name']}: {species['samples']} samples ({species['accession']})")
    else:
        print(f"  ⚠️  {species['name']}: Config not found")

# Save high-priority list
with open('output/ant_discovery/high_priority_species.json', 'w') as f:
    json.dump(high_priority, f, indent=2)

PYEOF

echo ""
echo "✅ High-priority configurations deployed"
echo ""

# Step 5: Verify deployed configurations
echo "Step 5: Verifying deployed configurations..."
echo ""
ls -lh config/amalgkit/amalgkit_*.yaml | tail -10
echo ""

# Step 6: Show example workflow commands
echo "Step 6: Example workflow commands..."
echo ""
echo "To run a single species workflow:"
echo ""
echo "  bash scripts/rna/amalgkit/run_amalgkit.sh \\"
echo "      --config config/amalgkit/amalgkit_camponotus_floridanus.yaml \\"
echo "      --steps metadata,select,getfastq,quant,merge,curate,sanity"
echo ""
echo "To run all high-priority species in parallel:"
echo ""
echo "  python3 scripts/rna/orchestrate_workflows.py \\"
echo "      --species camponotus_floridanus solenopsis_invicta \\"
echo "      --steps metadata select getfastq quant merge curate sanity"
echo ""
echo "To monitor progress:"
echo ""
echo "  python3 scripts/rna/get_current_status.py"
echo ""

# Step 7: Optional - Start workflows
echo "================================================================================"
echo "DISCOVERY AND DEPLOYMENT COMPLETE"
echo "================================================================================"
echo ""
echo "Next steps:"
echo "  1. Review: cat output/ant_discovery/DISCOVERY_REPORT.md"
echo "  2. Verify: ls -lh config/amalgkit/"
echo "  3. Execute: bash scripts/rna/amalgkit/run_amalgkit.sh --config ..."
echo ""
echo "To start workflows now, run:"
echo ""
echo "  bash scripts/rna/examples/run_high_priority_species.sh"
echo ""
echo "================================================================================"




