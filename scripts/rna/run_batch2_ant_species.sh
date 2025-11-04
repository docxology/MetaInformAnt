#!/usr/bin/env bash
# Run BATCH 2: Remaining 10 ant species (after batch 1 completes)

echo "================================================================================"
echo "LAUNCHING BATCH 2: REMAINING 10 ANT SPECIES"
echo "================================================================================"
echo ""
echo "Run this after Batch 1 (top 10) completes"
echo ""

# Batch 2: Remaining 10 species (total: 728 samples)
SPECIES=(
    "linepithema_humile"          # 173 samples
    "temnothorax_nylanderi"       # 124 samples
    "pogonomyrmex_barbatus"       # 120 samples
    "temnothorax_curvispinosus"   # 116 samples
    "acromyrmex_echinatior"       # 66 samples
    "lasius_neglectus"            # 51 samples
    "wasmannia_auropunctata"      # 33 samples
    "formica_exsecta"             # 23 samples
    "vollenhovia_emeryi"          # 15 samples
    "myrmica_rubra"               # 7 samples
)

echo "Species to process (Batch 2):"
echo " 11. Linepithema humile         - 173 samples"
echo " 12. Temnothorax nylanderi      - 124 samples"
echo " 13. Pogonomyrmex barbatus      - 120 samples"
echo " 14. Temnothorax curvispinosus  - 116 samples"
echo " 15. Acromyrmex echinatior      - 66 samples"
echo " 16. Lasius neglectus           - 51 samples"
echo " 17. Wasmannia auropunctata     - 33 samples"
echo " 18. Formica exsecta            - 23 samples"
echo " 19. Vollenhovia emeryi         - 15 samples"
echo " 20. Myrmica rubra              - 7 samples"
echo ""
echo "Total: 728 samples"
echo ""
echo "================================================================================"
echo ""

STARTED=0
FAILED=0

for species in "${SPECIES[@]}"; do
    config_file="config/amalgkit/amalgkit_${species}.yaml"
    
    if [ ! -f "$config_file" ]; then
        echo "⚠️  Config not found: $config_file"
        ((FAILED++))
        continue
    fi
    
    log_file="output/batch2_${species}_$(date +%Y%m%d_%H%M%S).log"
    
    echo "▶️  Starting: $species"
    echo "   Config: $config_file"
    echo "   Log: $log_file"
    
    nohup bash scripts/rna/amalgkit/run_amalgkit.sh \
        --config "$config_file" \
        --steps metadata,select,getfastq,quant,merge,curate,sanity \
        > "$log_file" 2>&1 &
    
    ((STARTED++))
    sleep 3
done

echo ""
echo "================================================================================"
echo "BATCH 2 LAUNCH COMPLETE"
echo "================================================================================"
echo ""
echo "✅ Started: $STARTED workflows"
if [ $FAILED -gt 0 ]; then
    echo "⚠️  Failed: $FAILED configs not found"
fi
echo ""
echo "Monitor progress:"
echo "  python3 scripts/rna/get_current_status.py"
echo "  tail -f output/batch2_*.log"
echo ""
echo "Estimated completion: 12-24 hours for batch 2"
echo "================================================================================"

