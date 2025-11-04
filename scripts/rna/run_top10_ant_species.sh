#!/usr/bin/env bash
# Run amalgkit workflows for TOP 10 ant species by sample count
# Focus: Complete end-to-end processing before moving to next batch

echo "================================================================================"
echo "LAUNCHING TOP 10 ANT SPECIES WORKFLOWS"
echo "================================================================================"
echo ""
echo "Strategy: Complete download + quant + merge + cleanup for 10 species"
echo "        Auto-delete FASTQs after quantification to save space"
echo ""

# Top 10 species by sample count (total: 3,820 samples)
SPECIES=(
    "harpegnathos_saltator"           # 695 samples
    "temnothorax_longispinosus"       # 557 samples
    "solenopsis_invicta"              # 451 samples
    "monomorium_pharaonis"            # 370 samples
    "camponotus_floridanus"           # 359 samples
    "temnothorax_rugatulus"           # 316 samples
    "ooceraea_biroi"                  # 278 samples
    "atta_cephalotes"                 # 239 samples
    "cardiocondyla_obscurior"         # 191 samples
    "lasius_niger"                    # 191 samples
)

echo "Species to process:"
echo "  1. Harpegnathos saltator      - 695 samples"
echo "  2. Temnothorax longispinosus  - 557 samples"
echo "  3. Solenopsis invicta         - 451 samples"
echo "  4. Monomorium pharaonis       - 370 samples"
echo "  5. Camponotus floridanus      - 359 samples"
echo "  6. Temnothorax rugatulus      - 316 samples"
echo "  7. Ooceraea biroi             - 278 samples"
echo "  8. Atta cephalotes            - 239 samples"
echo "  9. Cardiocondyla obscurior    - 191 samples"
echo " 10. Lasius niger               - 191 samples"
echo ""
echo "Total: 3,820 samples"
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
    
    log_file="output/top10_${species}_$(date +%Y%m%d_%H%M%S).log"
    
    echo "▶️  Starting: $species"
    echo "   Config: $config_file"
    echo "   Log: $log_file"
    
    # Run in background with nohup
    # Key: keep_fastq: no in configs ensures auto-cleanup
    nohup bash scripts/rna/amalgkit/run_amalgkit.sh \
        --config "$config_file" \
        --steps metadata,select,getfastq,quant,merge,curate,sanity \
        > "$log_file" 2>&1 &
    
    ((STARTED++))
    
    # Stagger starts to avoid overwhelming the system
    sleep 3
done

echo ""
echo "================================================================================"
echo "WORKFLOW LAUNCH COMPLETE"
echo "================================================================================"
echo ""
echo "✅ Started: $STARTED workflows"
if [ $FAILED -gt 0 ]; then
    echo "⚠️  Failed: $FAILED configs not found"
fi
echo ""
echo "Active workflows:"
sleep 2
ps aux | grep "run_amalgkit" | grep -v grep | wc -l
echo ""
echo "These workflows will:"
echo "  1. Download FASTQs for each sample"
echo "  2. Quantify with Kallisto (12 threads)"
echo "  3. Auto-delete FASTQs after quantification"
echo "  4. Merge expression matrices"
echo "  5. Run QC and validation"
echo ""
echo "Disk space management: FASTQs deleted immediately after quant"
echo ""
echo "Monitor progress:"
echo "  python3 scripts/rna/get_current_status.py"
echo "  tail -f output/top10_*.log"
echo ""
echo "Estimated completion: 24-48 hours for all 10 species"
echo "================================================================================"

