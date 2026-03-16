# 📊 MetaInformAnt Pipeline Progress Report

**Telemetry Refreshed:** 2026-03-16 03:30 UTC
**Orchestration Node:** `metainformant-pipeline` (n2-standard-16)
**Local Workstation Time:** 2026-03-15 20:30 local

---

## 📈 1. Live Progress Matrix (All Species)

| Species | Status | Quantified | Pending | Failed | Total | Downstream |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Apis mellifera** | 🔵 Running | 4,658 | 2,588 | 100 | 7,370 | ❌ Not run |
| **Ooceraea biroi** | ✅ Complete | 274 | 0 | 0 | 274 | ✅ Success |
| **Linepithema humile** | ✅ Complete | 173 | 0 | 0 | 173 | ✅ Success |
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
| **Wasmannia auropunctata** | ❌ Blocked | 33 | 0 | 0 | 33 | ❌ Curation Fail* |
| **Anoplolepis gracilipes** | ❌ Blocked | 7 | 0 | 0 | 7 | ❌ Curation Fail* |
| **Dinoponera quadriceps** | ❌ Blocked | 13 | 0 | 0 | 13 | ❌ Curation Fail* |
| **Vollenhovia emeryi** | ❌ Blocked | 15 | 0 | 0 | 15 | ❌ Curation Fail* |

*\*Downstream failures for small species are primarily due to "Metadata file not found" or "No samples meet filtering criteria" during the asynchronous curation loop. These will require manual reconciliation of LITE files.*

---

## 🧵 2. Thread Pool & Hardware Activity

The VM is an `n2-standard-16` with 16 vCPUs (hyperthreaded to 32 logical cores). We are operating with a **24-worker** thread pool for quantification and a dedicated core-set for downstream R analysis.

### Active Process Snapshot
| PID | Command | Resource | Activity |
| :--- | :--- | :--- | :--- |
| **403330** | `kallisto quant` | cpu: 63.4% | Quantifying `SRR25008587` (Apis) |
| **403867** | `kallisto quant` | cpu: 65.5% | Quantifying `SRR20272017` (Apis) |
| **397441** | `bg_curate_phase2.sh` | cpu: 0.1% | Idle (Phase 2 Loop Terminal) |
| **403150** | `pigz -p 2` | cpu: 104% | Decompressing `SRR...` buffer |

**Overall Saturation:** 🟢 **High Stability**
Quantification velocity is holding steady at **~123 samples/hour**. Currently at `4,658 / 7,370`.

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

## 🔭 4. Performance Assessment (T+03:30 UTC | 20:30 Local)

- **Curation Threshold**: 15 species are now fully finalized (`✅ Complete`). The Phase 2 curation successfully cleared `Ooceraea biroi` and `Linepithema humile`.
- **Apis Velocity**: Strong holding at **~123 samples/hour**. Currently at **4,658 / 7,370** finalized (63%).
- **Metadata Blocks**: Small species (`Wasmannia`, `Anoplolepis`, `Dinoponera`, `Vollenhovia`) are currently blocked by missing metadata files in the automated loop, likely requiring manual layout checks due to their small sample sizes.
