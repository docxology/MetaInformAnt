# 📊 MetaInformAnt Pipeline Progress Report

**Telemetry Refreshed:** 2026-03-16 12:15 UTC (T+~36h)
**Orchestration Node:** `metainformant-pipeline` (n2-standard-16)
**Local Workstation Time:** 2026-03-16 05:15 local

> [!IMPORTANT]
> **Pipeline Status: 🟢 HIGHSPEED QUANTIFICATION**
> `apis_mellifera` has crossed the **5,800** sample mark. Finalizing the remaining 1,400 samples at a velocity of **~128 samples/hour**.

---

## 📈 1. Live Progress Matrix (All Species)

| Species | Status | Quantified | Pending | Failed | Total | Downstream |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Apis mellifera** | 🔵 Running | 5,810 | 1,418 | 118 | 7,370 | ❌ Not run |
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
| PID | Command | Resource | Elapse (MM:SS) | Status |
| :--- | :--- | :--- | :--- | :--- |
| **838620** | `kallisto quant` | cpu: 60% | 47:56 | 🔵 Active |
| **841673** | `kallisto quant` | cpu: 65% | 26:48 | 🔵 Active |
| **843691** | `kallisto quant` | cpu: 62% | 10:24 | 🔵 Active |
| **844991** | `kallisto quant` | cpu: 68% | 00:29 | 🔵 Active |

**Performance Metrics:**
- **Quantification Velocity**: ~128 samples/hour (Current: `apis_mellifera`)
- **Disk IO Capacity**: 737GB available on master SSD.
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

## 🔭 4. Performance Assessment (T+~36h UTC)

- **Apis Achievement**: We are at **79%** of the final goal (**5,810 / 7,370**). Velocity remains ultra-stable at ~128/hr.
- **Curation Milestone**: 15 species are fully verified and finalized. Batch 2 added `Ooceraea biroi` and `Linepithema humile` to the Success list.
- **Targeting Completion**: At current velocity, the final 1,400 samples should clear in approximately ~11 hours.
