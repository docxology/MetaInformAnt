#!/usr/bin/env bash
# Run amalgkit workflows for all discovered ant species

echo "================================================================================"
echo "LAUNCHING AMALGKIT WORKFLOWS FOR ALL ANT SPECIES"
echo "================================================================================"
echo ""

# All newly discovered species with validated genomes (20 species)
SPECIES=(
    "harpegnathos_saltator"
    "temnothorax_longispinosus"
    "solenopsis_invicta"
    "monomorium_pharaonis"
    "camponotus_floridanus"
    "temnothorax_rugatulus"
    "ooceraea_biroi"
    "atta_cephalotes"
    "cardiocondyla_obscurior"
    "lasius_niger"
    "linepithema_humile"
    "temnothorax_nylanderi"
    "pogonomyrmex_barbatus"
    "temnothorax_curvispinosus"
    "acromyrmex_echinatior"
    "lasius_neglectus"
    "wasmannia_auropunctata"
    "formica_exsecta"
    "vollenhovia_emeryi"
    "myrmica_rubra"
)

echo "Starting workflows for ${#SPECIES[@]} ant species..."
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
    
    log_file="output/amalgkit_${species}_$(date +%Y%m%d_%H%M%S).log"
    
    echo "▶️  Starting: $species"
    echo "   Config: $config_file"
    echo "   Log: $log_file"
    
    # Run in background with nohup
    nohup bash scripts/rna/amalgkit/run_amalgkit.sh \
        --config "$config_file" \
        --steps metadata,select,getfastq,quant,merge,curate,sanity \
        > "$log_file" 2>&1 &
    
    ((STARTED++))
    
    # Small delay to stagger starts
    sleep 2
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
ps aux | grep "run_amalgkit" | grep -v grep | wc -l
echo ""
echo "Monitor progress:"
echo "  python3 scripts/rna/get_current_status.py"
echo ""
echo "Check active processes:"
echo "  ps aux | grep run_amalgkit | grep -v grep"
echo ""
echo "View logs:"
echo "  tail -f output/amalgkit_SPECIES_*.log"
echo ""
echo "================================================================================"

