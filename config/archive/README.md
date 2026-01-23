# Archived Configurations

This directory contains inactive or deprecated configuration files preserved for reference.

## Archived Species Configurations

These amalgkit configurations were moved from `config/amalgkit/` due to:
- Incomplete genome assemblies
- Low sample counts in SRA
- Pending validation
- Superseded by newer configurations

### Ant Species (Formicidae)

| Configuration | Species | Notes |
|---------------|---------|-------|
| `amalgkit_acromyrmex_echinatior.yaml` | *Acromyrmex echinatior* | Leafcutter ant |
| `amalgkit_atta_cephalotes.yaml` | *Atta cephalotes* | Leafcutter ant |
| `amalgkit_camponotus_floridanus.yaml` | *Camponotus floridanus* | Florida carpenter ant |
| `amalgkit_harpegnathos_saltator.yaml` | *Harpegnathos saltator* | Jerdon's jumping ant |
| `amalgkit_lasius_neglectus.yaml` | *Lasius neglectus* | Invasive garden ant |
| `amalgkit_monomorium_pharaonis.yaml` | *Monomorium pharaonis* | Pharaoh ant |
| `amalgkit_myrmica_rubra.yaml` | *Myrmica rubra* | European fire ant |
| `amalgkit_ooceraea_biroi.yaml` | *Ooceraea biroi* | Clonal raider ant |
| `amalgkit_solenopsis_invicta.yaml` | *Solenopsis invicta* | Red fire ant |
| `amalgkit_temnothorax_curvispinosus.yaml` | *Temnothorax curvispinosus* | Acorn ant |
| `amalgkit_temnothorax_longispinosus.yaml` | *Temnothorax longispinosus* | Acorn ant |
| `amalgkit_temnothorax_rugatulus.yaml` | *Temnothorax rugatulus* | Acorn ant |
| `amalgkit_vollenhovia_emeryi.yaml` | *Vollenhovia emeryi* | |
| `amalgkit_wasmannia_auropunctata.yaml` | *Wasmannia auropunctata* | Little fire ant |

### Other Species

| Configuration | Species | Notes |
|---------------|---------|-------|
| `amalgkit_amellifera.yaml` | *Apis mellifera* | Western honey bee |

## Reactivating Configurations

To reactivate an archived configuration:

1. Copy to `config/amalgkit/`:
   ```bash
   cp config/archive/amalgkit_<species>.yaml config/amalgkit/
   ```

2. Verify genome assembly is still current at NCBI

3. Update FTP URLs if needed

4. Test with:
   ```bash
   python3 scripts/rna/run_workflow.py --config config/amalgkit/<config>.yaml --check
   ```

## Notes

- Archived configs may reference outdated genome assemblies
- FTP URLs may have changed since archiving
- Sample availability in SRA may have changed
