#!/bin/bash
# Migrate legacy amellifera data to new structure

SPECIES="amellifera"
BASE="blue/amalgkit/$SPECIES"
OLD_QUANT="$BASE/quant"
NEW_WORK="$BASE/work"
NEW_QUANT="$NEW_WORK/quant"
NEW_FASTQ="$BASE/fastq/getfastq"

echo "Migrating $SPECIES data..."

# Ensure target directories exist
mkdir -p "$NEW_QUANT"
mkdir -p "$NEW_FASTQ"

# 1. Migrate Quant Results (Direct SRR folders in output root)
if [ -d "$OLD_QUANT" ]; then
    echo "Found old quant dir: $OLD_QUANT"
    
    # Move SRR directories directly under old quant
    # Find directories starting with SRR/ERR/DRR
    find "$OLD_QUANT" -maxdepth 1 -type d -name "[SED]RR*" -exec mv -v {} "$NEW_QUANT/" \;
    
    # 2. Migrate Nested Quant (quant/quant/SRR...)
    if [ -d "$OLD_QUANT/quant" ]; then
        echo "Found nested quant dir: $OLD_QUANT/quant"
        find "$OLD_QUANT/quant" -maxdepth 1 -type d -name "[SED]RR*" -exec mv -v {} "$NEW_QUANT/" \;
    fi
    
    # 3. Migrate Fastq (quant/getfastq/...)
    # Sometimes it's getfastq directly or inside quant
    if [ -d "$OLD_QUANT/getfastq" ]; then
        echo "Found old fastq dir: $OLD_QUANT/getfastq"
        find "$OLD_QUANT/getfastq" -maxdepth 1 -type d -name "[SED]RR*" -exec mv -v {} "$NEW_FASTQ/" \;
    fi
else
    echo "No legacy quant directory found at $OLD_QUANT"
fi

echo "Migration complete."
echo "New Quant Count: $(find "$NEW_QUANT" -maxdepth 1 -type d | wc -l)"
