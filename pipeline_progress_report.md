# 🧬 MetaInformAnt Orchestration: GCP Telemetry Dashboard

This document provides a highly visual, accessible framework for monitoring the end-to-end `amalgkit` RNA-seq pipeline on the GCP compute node. It covers both the current automated progress metrics and precisely **how** you can reproduce these telemetry pulls effortlessly from your local terminal.

> [!TIP]
> **Zero-Interruption Monitoring Policy**
> The commands documented below natively interrogate the active container orchestrator WITHOUT stopping or slowing down the primary ThreadPool. 

---

## 📈 1. Live Progress Matrix (All Species)

*Last Polled: 2026-03-15 15:40 UTC*

### The Metric Command
Safely dump the master SQLite pipeline matrix using the bespoke orchestration script inside the running container.
```bash
gcloud compute ssh metainformant-pipeline \
    --zone=us-central1-a \
    --project=cryptoptera \
    --tunnel-through-iap \
    --ssh-flag="-q" \
    --command='sudo docker exec metainformant-pipeline-fresh python3 scripts/rna/check_pipeline_status.py -v'
```

### Current Status
```text
================================================================================
  Amalgkit Pipeline Status  │  Process: 🟢 RUNNING
================================================================================
Species                       pending  downloa  downloa  quantif  quantif   failed   Total  Downstream      
--------------------------------------------------------------------------------
  anoplolepis_gracilipes            0        0        0        0        7        0       7  ❌ Not run
  acromyrmex_echinatior             0        0        0        0       44        0      44  ✅ Complete
  dinoponera_quadriceps             0        0        0        0       13        0      13  ❌ Not run
  vollenhovia_emeryi                0        0        0        0       15        0      15  ❌ Not run
  odontomachus_brunneus             0        0        0        0       19        0      19  ✅ Complete
  formica_exsecta                   0        0        0        0       23        0      23  ⚠️  Merge only
  temnothorax_americanus            0        0        0        0      331        0     331  ❌ Not run
  wasmannia_auropunctata            0        0        0        0       33        0      33  ❌ Not run
  nylanderia_fulva                  0        0        0        0       40        0      40  ⚠️  Merge only
  temnothorax_curvispinosus         0        0        0        0       43        0      43  ⚠️  Merge only
  pogonomyrmex_barbatus             0        0        0        0      132        0     132  ❌ Not run
  cardiocondyla_obscurior           0        0        0        0      167        0     167  ⚠️  Merge only
  temnothorax_nylanderi             0        0        0        0      154       12     166  ⚠️  Merge only
  linepithema_humile                0        0        0        0      173        0     173  ❌ Not run
  atta_cephalotes                   0        0        0        0      217        0     217  ✅ Complete
  ooceraea_biroi                    0        0        0        0      274        0     274  ❌ Not run
  camponotus_floridanus             0        0        0        0      366        1     367  ⚠️  Merge only
  solenopsis_invicta                0        0        0        0      450        1     451  ✅ Complete
  monomorium_pharaonis              0        0        0        0      370        0     370  ⚠️  Merge only
  temnothorax_longispinosus         0        0        0        0      508        0     508  ✅ Complete
  harpegnathos_saltator             0        0        0        0      689        0     689  ✅ Complete
  apis_mellifera                 4132        3        0       21     3143       71    7370  ❌ Not run
--------------------------------------------------------------------------------
  TOTAL                          4132        3        0       21     7211       85   11452
================================================================================
```
================================================================================
```
### 🐝 Apis Mellifera Velocity Log
- **T+~65 min** (00:50 UTC Day 3): `1,331 / 7,370` Complete (+178 burst)
- **T+~4 hours** (05:07 UTC Day 3): `1,879 / 7,370` Complete (+548 huge files flawlessly processed)
- **T+~11.5 hours** (12:20 UTC Day 3): `2,817 / 7,370` Complete (+938 natively processed overnight by the 24-core local SSD engine)
- **T+~12 hours** (12:44 UTC Day 3): `2,877 / 7,370` Complete (+60 heavy SRAs quantified natively trailing 24 minutes)
- **T+~12.5 hours** (13:10 UTC Day 3): `2,899 / 7,370` Complete (+22 quantified natively trailing 26 minutes)
- **T+~14.5 hours** (15:08 UTC Day 3): `3,098 / 7,370` Complete (+199 quantified natively trailing 2 hours)
- **T+~15 hours** (15:40 UTC Day 3): `3,143 / 7,370` Complete (+45 quantified natively trailing 32 minutes)

The ThreadPool is maintaining a consistent velocity of **~85-100 samples/hour**. Since the 13:10 UTC poll, we have added **+244 processed samples** to the `apis_mellifera` count.

---

### 🔍 Current Activity Assessment (T+12:30 local / 15:40 UTC)
- **Quantification Processes (`quant`)**: All 24 cores remain fully saturated. Velocity is averaging ~90 samples/hour as the ThreadPool processes `apis_mellifera`.
- **Post-Merge Processes (`downstream`)**: The `sva` dependency patch successfully unblocked the curation engine. **3 species** (`odontomachus_brunneus`, `solenopsis_invicta`, and `temnothorax_longispinosus`) have already moved from `⚠️ Merge only` to `✅ Complete`. The background loop is currently processing the remaining 7 species in the backlog.

---

## ⚙️ 2. Hardware Resource & Thread Pool Tracking

Confirm that CPU saturation and bandwidth are fully utilized by examining the real background processes. This bypasses the shell to show exactly what `multiprocessing` is doing.

### The Diagnostic Command
```bash
gcloud compute ssh metainformant-pipeline \
    --zone=us-central1-a \
    --project=cryptoptera \
    --tunnel-through-iap \
    --ssh-flag="-q" \
    --command='echo -e "\n--- THREAD POOL ---\n" && sudo docker top metainformant-pipeline-fresh'
```

### Current Status (Docker ps)
```text
UID                 PID                 PPID                C                   STIME               TTY                 TIME                CMD
root                245176              245157              0                   Mar13               ?                   00:00:00            /bin/bash /app/run_pipeline.sh
root                245218              245176              0                   Mar13               ?                   00:01:46            python3 scripts/rna/run_all_species.py --max-gb 350.0 --workers 24 --threads 2

< - - - 24 Independent Simultaneous Thread Pools of NCBI & Kallisto - - - >

root                622404              245218              67.0                12:43               ?                   17:40               kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR13871706 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR13871706/SRR13871706_1.fastq.gz
root                622650              245218              66.7                12:45               ?                   16:43               kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR37117705 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR37117705/SRR37117705_1.fastq.gz
root                622910              245218              62.5                12:47               ?                   14:25               kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR25389649 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR25389649/SRR25389649_1.fastq.gz
root                622958              245218              65.5                12:47               ?                   14:51               kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR29535100 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29535100/SRR29535100_1.fastq.gz

< - - - Background Single-Thread Curation Engine - - - >

root                296743              1                   0                   12:12               ?                   00:00:00            /bin/bash scripts/rna/bg_curate.sh
root                296842              296744              5                   12:13               ?                   00:00:00            /app/.venv/bin/python3 /app/.venv/bin/amalgkit merge --out_dir output/amalgkit/temnothorax_longispinosus/merged --metadata output/amalgkit/temnothorax_longispinosus/work/metadata/metadata_selected.tsv
```

---

## 🛠️ 3. Spawning Detached Asynchronous Curation (`Merge only` Backlog)

The main 24 cores are entirely dedicated to `apis_mellifera` quantification logic in Phase 1. 

To clear the backlog of `⚠️ Merge only` species cleanly without halting quantification bandwidth, we inject an independent single-thread bash script into the orchestration image.

### The Injection Sequence
```bash
# 1. Create a looping bash script locally
cat << 'EOF' > bg_curate.sh
#!/bin/bash
for sp in temnothorax_longispinosus solenopsis_invicta monomorium_pharaonis camponotus_floridanus cardiocondyla_obscurior temnothorax_nylanderi temnothorax_curvispinosus nylanderia_fulva formica_exsecta odontomachus_brunneus; do 
    echo "Processing $sp"
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_${sp}.yaml --no-progress --steps merge curate sanity
done
EOF

# 2. Transfer script to the VM host natively via IAP
gcloud compute scp bg_curate.sh metainformant-pipeline:~/ --zone=us-central1-a --project=cryptoptera --tunnel-through-iap

# 3. Inject and run detached in the orchestrator container
gcloud compute ssh metainformant-pipeline \
    --zone=us-central1-a \
    --project=cryptoptera \
    --tunnel-through-iap \
    --ssh-flag="-q" \
    --command='sudo docker cp ~/bg_curate.sh metainformant-pipeline-fresh:/app/scripts/rna/bg_curate.sh && sudo docker exec metainformant-pipeline-fresh bash -c "chmod +x scripts/rna/bg_curate.sh" && sudo docker exec -d metainformant-pipeline-fresh bash -c "nohup scripts/rna/bg_curate.sh > output/amalgkit/manual_downstream.log 2>&1 &"'
```

### Checking the Asynchronous Curation Log
Monitor exactly where the curation engine is failing, passing, or processing at any time.

```bash
gcloud compute ssh metainformant-pipeline \
    --zone=us-central1-a \
    --project=cryptoptera \
    --tunnel-through-iap \
    --ssh-flag="-q" \
    --command='sudo docker exec metainformant-pipeline-fresh tail -n 20 output/amalgkit/manual_downstream.log'
```
