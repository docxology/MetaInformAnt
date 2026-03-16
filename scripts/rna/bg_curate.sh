#!/bin/bash
for sp in temnothorax_longispinosus solenopsis_invicta monomorium_pharaonis camponotus_floridanus cardiocondyla_obscurior temnothorax_nylanderi temnothorax_curvispinosus nylanderia_fulva formica_exsecta odontomachus_brunneus; do 
    echo "Processing $sp"
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_${sp}.yaml --no-progress --steps merge curate sanity
done
