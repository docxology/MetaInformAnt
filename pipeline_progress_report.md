# 📊 MetaInformAnt Pipeline Progress Report — Recovery VM

**Telemetry Refreshed:** 2026-03-22 00:30 UTC  
**Orchestration Node:** `metainformant-recovery` (n2-standard-32, us-central1-a)  
**Local Workstation Time:** 2026-03-21 17:30 local

> [!NOTE]
> **Pipeline Status: 🟢 RECOVERING — polistes_fuscatus actively quantifying**  
> Emergency disk cleanup freed **213+ GB**. Disk now ~63 GB free. All failed samples re-queued to pending.

---

## 🛠️ Fixes Applied This Session

| Fix | Status | Detail |
| :--- | :---: | :--- |
| FASTQ cleanup (post-quant) | ✅ Done | 213+ GB freed (141+ files deleted, daemon running) |
| Failed samples re-queued | ✅ Done | 139 samples reset to pending across two rounds |
| Disk throttle patch (source) | ✅ Done | `streaming_orchestrator.py`: 10 GB → 50 GB on VM |
| Cleanup daemon v2 launched | ✅ Done | Runs every 30s, cleans completed FASTQs + re-queues failures |

---

## 🌿 Wave-2 Species Progress (Recovery VM)

| Species | Quantifying | Downloading | Pending | Failed | Total | % Done |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Polistes canadensis** | 0 | 0 | 115 | 0 | 115 | ~87% (100 on disk) |
| **Polistes fuscatus** | 19 | 1 | ~260 | ~0* | 293 | ~12% |
| **Megachile rotundata** | 0 | 0 | 170 | 0 | 170 | Not started |
| **Athalia rosae** | 0 | 0 | 43 | 0 | 43 | Not started |
| **Apis mellifera** (recovery) | 0 | 0 | 7,264 | 0 | 7,264 | Not started |

*\*Failures continuously re-queued by daemon every 30s*

> **Disk**: ~63 GB free (was 7.3 GB). Daemon v2 cleaning FASTQs every 30s.  
> **Workers**: 20 workers x 2 threads = 40 concurrent threads.

---

## 🏆 Wave-1 Species — Previously Completed

| Species | Quantified | Failed | Total |
| :--- | :---: | :---: | :---: |
| Apis mellifera | 7,222 | 148 | 7,370 |
| Harpegnathos saltator | 689 | 0 | 689 |
| Temnothorax longispinosus | 508 | 0 | 508 |
| Solenopsis invicta | 450 | 1 | 451 |
| Monomorium pharaonis | 370 | 0 | 370 |
| Camponotus floridanus | 366 | 1 | 367 |
| Temnothorax americanus | 331 | 0 | 331 |
| Ooceraea biroi | 274 | 0 | 274 |
| Atta cephalotes | 217 | 0 | 217 |
| Linepithema humile | 173 | 0 | 173 |
| Cardiocondyla obscurior | 167 | 0 | 167 |
| Temnothorax nylanderi | 154 | 12 | 166 |
| Pogonomyrmex barbatus | 132 | 0 | 132 |
| *(+ 8 smaller species)* | 233 | 0 | 233 |
| **TOTAL Wave-1** | **11,278** | **162** | **11,440** |

---

## 🔧 Monitoring Commands

```bash
# Quick DB state check
gcloud compute ssh metainformant-recovery --zone=us-central1-a \
  --command='sqlite3 -readonly /opt/MetaInformAnt/output/amalgkit/pipeline_progress.db \
  "SELECT species,state,COUNT(*) FROM samples GROUP BY species,state ORDER BY species,state;"'

# Daemon v2 log (cleanup + re-queue status)
gcloud compute ssh metainformant-recovery --zone=us-central1-a \
  --command='tail -20 /tmp/daemon2.log'

# Disk free
gcloud compute ssh metainformant-recovery --zone=us-central1-a \
  --command='df -h / | tail -1'

# Re-queue any new failures manually
gcloud compute ssh metainformant-recovery --zone=us-central1-a \
  --command='sudo sqlite3 /opt/MetaInformAnt/output/amalgkit/pipeline_progress.db \
  "UPDATE samples SET state=\"pending\",error=NULL WHERE state=\"failed\"; SELECT changes();"'
```

---

## ⚡ Root Cause & Fix Summary

| Problem | Root Cause | Fix |
| :--- | :--- | :--- |
| Disk 100% full | FASTQ files not deleted after quantification | Cleanup daemon (running every 30s) |
| `raise RuntimeError()` failures | `curl: (23)` write error — disk full during download | FASTQ cleanup + throttle raised to 50GB |
| Workers all throttle-paused | Threshold 10 GB too low for 20 concurrent downloaders | Patched source to 50 GB on VM |
| Failed samples stuck | Error state not auto-cleared | Daemon v2 re-queues failures every 30s |
