# 📊 MetaInformAnt Pipeline Progress Report

**Telemetry Refreshed:** 2026-03-16 17:48 UTC (T+~41.5h)
**Orchestration Node:** `metainformant-pipeline` (n2-standard-16)
**Local Workstation Time:** 2026-03-16 10:48 local

> [!IMPORTANT]
> **Pipeline Status: 🟢 HIGHSPEED QUANTIFICATION**
> `apis_mellifera` has crossed the **6,400** sample mark. Finalizing the remaining 744 samples at a velocity of **~120-125 samples/hour**.

---

## 📈 1. Live Progress Matrix (All Species)

| Species | Status | Quantified | Pending | Failed | Total | Downstream |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Apis mellifera** | 🔵 Running | 6,468 | 744 | 134 | 7,370 | ❌ Not run |
| **Ooceraea biroi** | ✅ Complete | 274 | 0 | 0 | 274 | ✅ Success |
| **Linepithema humile** | ✅ Complete | 173 | 0 | 0 | 173 | ✅ Success |
| **Pogonomyrmex barbatus** | ✅ Complete | 132 | 0 | 0 | 132 | ✅ Success |
| **Temnothorax americanus** | ✅ Complete | 331 | 0 | 0 | 331 | ✅ Success |
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
| **Wasmannia auropunctata** | ❌ Blocked | 33 | 0 | 0 | 33 | ❌ Curation Fail* |
| **Anoplolepis gracilipes** | ❌ Blocked | 7 | 0 | 0 | 7 | ❌ Curation Fail* |
| **Dinoponera quadriceps** | ❌ Blocked | 13 | 0 | 0 | 13 | ❌ Curation Fail* |
| **Vollenhovia emeryi** | ❌ Blocked | 15 | 0 | 0 | 15 | ❌ Curation Fail* |

*\*Downstream failures for the smallest species are due to metadata filtering constraints (e.g., LITE-only runs). These will be resolved during final post-pipeline cleanup.*

---

## 🧵 2. Thread Pool & Hardware Saturation

The VM is operating under a **24-worker** quantification pool. CPU saturation is sustainable, and SSD buffer space is optimal.

### Active Process Snapshot (Live Runtimes)
| PID | Command | Resource | Elapse (HH:MM:SS) | Status |
| :--- | :--- | :--- | :--- | :--- |
| **501081** | `kallisto quant` | cpu: 100% | 04:33:21 | 🔵 Active (Heavy) |
| **524494** | `kallisto quant` | cpu: 100% | 00:28:17 | 🔵 Active |
| **525924** | `kallisto quant` | cpu: 100% | 00:08:19 | 🔵 Active |
| **526435** | `kallisto quant` | cpu: 100% | 00:05:20 | 🔵 Active |

**Performance Metrics:**
- **Quantification Velocity**: ~120 samples/hour (Current: `apis_mellifera`)
- **Disk IO Capacity**: Stable on master SSD.
- **Background Curation**: Phase 2 loop is completed.

---

## 🛠️ 3. Diagnostic Command Toolkit

Run these commands inside your terminal to fetch the latest live data from the GCP orchestrator.

### Fetch Status Matrix
```bash
docker exec metainformant-pipeline-fresh python3 scripts/rna/check_pipeline_status.py -v
```

### Audit Active Compute
```bash
docker exec metainformant-pipeline-fresh ps -eo pid,cmd,etime,stat | grep kallisto | grep -v grep
```

### Tail Phase 2 Logs
```bash
docker exec metainformant-pipeline-fresh tail -n 50 output/amalgkit/manual_downstream_phase2.log
```

---

## 🔭 4. Performance Assessment (T+~41.5h UTC)

- **Apis Achievement**: We are at **87.7%** of the final goal (**6,468 / 7,370**). Velocity remains robust despite some heavier sample outliers.
- **Curation Milestone**: 15 species are fully verified and finalized. Batch 2 added `Ooceraea biroi` and `Linepithema humile` to the Success list.
- **Targeting Completion**: At current velocity, the final 744 samples should clear in approximately **~6 hours**.
