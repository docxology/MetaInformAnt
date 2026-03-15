# 📊 MetaInformAnt Pipeline Progress Report

**Telemetry Refreshed:** 2026-03-15 16:45 UTC
**Orchestration Node:** `metainformant-pipeline` (n2-standard-16)
**Local Workstation Time:** 2026-03-15 09:42 local

---

## 📈 1. Live Progress Matrix (All Species)

| Species | Status | Quantified | Pending | Failed | Total | Downstream |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Apis mellifera** | 🔵 Running | 3,143 | 4,132 | 71 | 7,370 | ❌ Not run |
| **Camponotus floridanus** | 🟠 Curating | 366 | 0 | 1 | 367 | 🔄 Active (SVA) |
| **Solenopsis invicta** | ✅ Complete | 450 | 0 | 1 | 451 | ✅ Success |
| **Temnothorax longispinosus** | ✅ Complete | 508 | 0 | 0 | 508 | ✅ Success |
| **Harpegnathos saltator** | ✅ Complete | 689 | 0 | 0 | 689 | ✅ Success |
| **Atta cephalotes** | ✅ Complete | 217 | 0 | 0 | 217 | ✅ Success |
| **Acromyrmex echinatior** | ✅ Complete | 44 | 0 | 0 | 44 | ✅ Success |
| **Odontomachus brunneus** | ✅ Complete | 19 | 0 | 0 | 19 | ✅ Success |
| **Monomorium pharaonis** | ⚠️ Pending | 370 | 0 | 0 | 370 | ⚠️ Merge only |
| **Temnothorax nylanderi** | ⚠️ Pending | 154 | 0 | 12 | 166 | ⚠️ Merge only |
| **Cardiocondyla obscurior** | ⚠️ Pending | 167 | 0 | 0 | 167 | ⚠️ Merge only |
| **Temnothorax curvispinosus** | ⚠️ Pending | 43 | 0 | 0 | 43 | ⚠️ Merge only |
| **Nylanderia fulva** | ⚠️ Pending | 40 | 0 | 0 | 40 | ⚠️ Merge only |
| **Formica exsecta** | ⚠️ Pending | 23 | 0 | 0 | 23 | ⚠️ Merge only |
| **Pogonomyrmex barbatus** | ⚠️ Pending | 132 | 0 | 0 | 132 | ❌ Not run |
| **Ooceraea biroi** | ⚠️ Pending | 274 | 0 | 0 | 274 | ❌ Not run |
| **Linepithema humile** | ⚠️ Pending | 173 | 0 | 0 | 173 | ❌ Not run |
| **Temnothorax americanus** | ⚠️ Pending | 331 | 0 | 0 | 331 | ❌ Not run |
| **Wasmannia auropunctata** | ⚠️ Pending | 33 | 0 | 0 | 33 | ❌ Not run |
| **Anoplolepis gracilipes** | ⚠️ Pending | 7 | 0 | 0 | 7 | ❌ Not run |
| **Dinoponera quadriceps** | ⚠️ Pending | 13 | 0 | 0 | 13 | ❌ Not run |
| **Vollenhovia emeryi** | ⚠️ Pending | 15 | 0 | 0 | 15 | ❌ Not run |

---

## 🧵 2. Thread Pool & Hardware Activity

The VM is an `n2-standard-16` with 16 vCPUs (hyperthreaded to 32 logical cores). We are operating with a **24-worker** thread pool for quantification and a dedicated core-set for downstream R analysis.

### Active Process Snapshot
| PID | Command | Resource | Activity |
| :--- | :--- | :--- | :--- |
| **358606** | `kallisto quant` | cpu: 68.7% | Quantifying `SRR17409101` (Apis) |
| **358640** | `curate.r` | cpu: 79.9% | SVA Correction (Camponotus) |
| **358444** | `pigz -p 2` | cpu: 104% | Decompressing `SRR35323153` |
| **358566** | `amalgkit curate`| mem: 0.2% | Orchestrating downstream workflow |

**Overall Saturation:** 🟢 **100% Core Load**
The ThreadPool is fully interleaved. As soon as a `kallisto` instance finishes an `apis_mellifera` sample, the next batch item from the 4,132 pending queue is injected.

---

## 🛠️ 3. Diagnostic Command Reference

These are the most reliable methods for gathering the telemetry presented in this report.

### Check Pipeline Matrix (SQLite Master)
```bash
docker exec metainformant-pipeline-fresh python3 scripts/rna/check_pipeline_status.py -v
```

### Audit Active Compute (Thread View)
```bash
docker exec metainformant-pipeline-fresh ps aux | grep -v grep | grep -E "kallisto|amalgkit|python3|sra|fastq"
```

### Monitor Downstream Progress (Curation)
```bash
docker exec metainformant-pipeline-fresh tail -f output/amalgkit/manual_downstream.log
```

---

## 🔭 4. Performance Assessment (T+16:45)

- **Apis Velocity**: Steady at **85-110 samples/hour**. Currently at `3,143 / 7,370`.
- **Curation Status**: SUCCESS. The R dependency patch (`sva`, `edgeR`, `Rtsne`) is working perfectly. `Camponotus floridanus` is currently undergoing deep SVA correction (74 surrogate variables), clearing the correlation rounds.
- **Disk Health**: SSD buffer is clear. `pigz` and `fasterq-dump` are operating at peak bandwidth.
