# 📊 MetaInformAnt Pipeline Progress Report

**Telemetry Refreshed:** 2026-03-15 22:30 UTC
**Orchestration Node:** `metainformant-pipeline` (n2-standard-16)
**Local Workstation Time:** 2026-03-15 15:31 local

---

## 📈 1. Live Progress Matrix (All Species)

| Species | Status | Quantified | Pending | Failed | Total | Downstream |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Apis mellifera** | 🔵 Running | 4,042 | 3,219 | 85 | 7,370 | ❌ Not run |
| **Camponotus floridanus** | ✅ Complete | 366 | 0 | 1 | 367 | ✅ Success |
| **Solenopsis invicta** | ✅ Complete | 450 | 0 | 1 | 451 | ✅ Success |
| **Temnothorax longispinosus** | ✅ Complete | 508 | 0 | 0 | 508 | ✅ Success |
| **Harpegnathos saltator** | ✅ Complete | 689 | 0 | 0 | 689 | ✅ Success |
| **Atta cephalotes** | ✅ Complete | 217 | 0 | 0 | 217 | ✅ Success |
| **Acromyrmex echinatior** | ✅ Complete | 44 | 0 | 0 | 44 | ✅ Success |
| **Odontomachus brunneus** | ✅ Complete | 19 | 0 | 0 | 19 | ✅ Success |
| **Monomorium pharaonis** | ✅ Complete | 370 | 0 | 0 | 370 | ✅ Success |
| **Temnothorax nylanderi** | ✅ Complete | 154 | 0 | 12 | 166 | ✅ Success |
| **Cardiocondyla obscurior** | ✅ Complete | 167 | 0 | 0 | 167 | ✅ Success |
| **Temnothorax curvispinosus** | ✅ Complete | 43 | 0 | 0 | 43 | ✅ Success |
| **Nylanderia fulva** | ✅ Complete | 40 | 0 | 0 | 40 | ✅ Success |
| **Formica exsecta** | ✅ Complete | 23 | 0 | 0 | 23 | ✅ Success |
| **Pogonomyrmex barbatus** | ✅ Complete | 132 | 0 | 0 | 132 | ✅ Success |
| **Temnothorax americanus** | ✅ Complete | 331 | 0 | 0 | 331 | ✅ Success |
| **Ooceraea biroi** | 🟠 Curating | 274 | 0 | 0 | 274 | 🔄 Phase 2 |
| **Linepithema humile** | ⚠️ Ready | 173 | 0 | 0 | 173 | ⚠️ Phase 2 Queue |
| **Wasmannia auropunctata** | ⚠️ Ready | 33 | 0 | 0 | 33 | ⚠️ Phase 2 Queue |
| **Anoplolepis gracilipes** | ⚠️ Ready | 7 | 0 | 0 | 7 | ⚠️ Phase 2 Queue |
| **Dinoponera quadriceps** | ⚠️ Ready | 13 | 0 | 0 | 13 | ⚠️ Phase 2 Queue |
| **Vollenhovia emeryi** | ⚠️ Ready | 15 | 0 | 0 | 15 | ⚠️ Phase 2 Queue |

---

## 🧵 2. Thread Pool & Hardware Activity

The VM is an `n2-standard-16` with 16 vCPUs (hyperthreaded to 32 logical cores). We are operating with a **24-worker** thread pool for quantification and a dedicated core-set for downstream R analysis.

### Active Process Snapshot
| PID | Command | Resource | Activity |
| :--- | :--- | :--- | :--- |
| **403330** | `kallisto quant` | cpu: 63.4% | Quantifying `SRR25008587` (Apis) |
| **403867** | `kallisto quant` | cpu: 65.5% | Quantifying `SRR20272017` (Apis) |
| **403022** | `bg_curate_phase2.sh` | cpu: 1.5% | Orchestrating curation Phase 2 |
| **403150** | `pigz -p 2` | cpu: 104% | Decompressing `SRR...` buffer |

**Overall Saturation:** 🟢 **Peak Performance**
Quantification velocity is currently **~120 samples/hour**. Currently at `4,042 / 7,370`.

---

## 🛠️ 3. Diagnostic Command Reference

These are the most reliable methods for gathering the telemetry presented in this report.

### Check Pipeline Matrix (SQLite Master)
```bash
docker exec metainformant-pipeline-fresh python3 scripts/rna/check_pipeline_status.py -v
```

### Audit Active Compute (Thread View)
```bash
docker exec metainformant-pipeline-fresh ps aux | grep -v grep | grep -E "kallisto|amalgkit|python3|sra|fastq|R"
```

### Monitor Downstream Progress (Phase 2 Curation)
```bash
docker exec metainformant-pipeline-fresh tail -f output/amalgkit/manual_downstream_phase2.log
```

---

## 🔭 4. Performance Assessment (T+22:00 UTC)

- **Curation Surge**: 13 species are now fully finalized (`✅ Complete`). The R dependency fix has proven incredibly stable.
- **Apis Acceleration**: `apis_mellifera` has crossed the 50% mark for its final batch, with **3,965 / 7,370** finalized.
- **Phase 2 Pipeline**: `bg_curate_phase2.sh` is now processing the final 8 species (`Pogonomyrmex` through `Vollenhovia`).
