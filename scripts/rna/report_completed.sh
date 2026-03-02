echo -e "| Species | Completed Samples | Total SRA Samples |"
echo -e "|---|---|---|"
for cfg in config/amalgkit/amalgkit_*.yaml; do
  species=$(basename $cfg .yaml | sed 's/amalgkit_//')
  if [[ "$species" == *"template"* || "$species" == *"test"* || "$species" == *"cross_species"* || "$species" == *"apis_mellifera_all"* ]]; then continue; fi
  
  total=0
  meta="output/amalgkit/$species/work/metadata/metadata.tsv"
  if [ -f "$meta" ]; then
    total=$(($(wc -l < $meta) - 1))
    if [ "$total" -lt 0 ]; then total=0; fi
  fi
  
  comp=$(find "output/amalgkit/$species" -name "*_abundance.tsv" 2>/dev/null | wc -l)
  
  if [ "$total" -gt 0 ] || [ "$comp" -gt 0 ]; then
    echo "| \`$species\` | **$comp** | $total |"
  fi
done
