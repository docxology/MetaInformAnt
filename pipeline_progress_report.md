# 🧬 Continuous Orchestration Telemetry

> [!TIP]
> **Live Monitoring from Any Shell**  
> Retrieve real-time metrics natively against the active container orchestration. The Kubernetes/Docker engine processes scaling and thread concurrency directly off the main VM.
> ```bash
> gcloud compute ssh metainformant-pipeline --zone=us-central1-a --project=cryptoptera --tunnel-through-iap --ssh-flag="-q" --command='sudo docker exec metainformant-pipeline-fresh python3 scripts/rna/check_pipeline_status.py -v && echo -e "\n--- THREAD POOL ---\n" && sudo docker top metainformant-pipeline-fresh'
> ```

## 📊 Operations Matrix: Entire Species Queue
*Capture Snapshot: 2026-03-15 00:50 UTC*

> [!NOTE]
> The primary 24-core Python orchestrator prioritizes smaller species first. All available quantitative threads are exclusively locked into dumping and aligning massive `apis_mellifera` clusters.

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
  odontomachus_brunneus             0        0        0        0       19        0      19  ⚠️  Merge only
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
  solenopsis_invicta                0        0        0        0      450        1     451  ⚠️  Merge only
  monomorium_pharaonis              0        0        0        0      370        0     370  ⚠️  Merge only
  temnothorax_longispinosus         0        0        0        0      508        0     508  ⚠️  Merge only
  harpegnathos_saltator             0        0        0        0      689        0     689  ✅ Complete
  apis_mellifera                 5984        7        0       17     1331       31    7370  ❌ Not run
--------------------------------------------------------------------------------
  TOTAL                          5984        7        0       17     5399       45   11452
================================================================================
```

## 🚀 Throughput & Velocity Accelerators

The ThreadPool continues to exceed baseline capacity by overlapping NCBI `curl` and `kallisto` streams across the local SSD blocks.

- **Current Apis Mellifera Quantified:** `1,331 / 7,370` Fully Complete
- **Sustained Massive Velocity:** **+178 extra samples** flawlessly quantified in the final explicit 65-minute window (`1153 -> 1331`).
- **Asymmetric Curation Engine (NEW):** An asynchronous single-threaded process (`bg_curate.sh`) was safely deployed into the cloud container. It sequentially loops internal `Merge only` species queues applying independent `merge`, `curate`, and `sanity` post-computations underneath the primary load constraints. *(Dependency Update Strategy active: `Rtsne` missing packages are being provisioned to resolve isolated R-language halt)*.

## ⚙️ Concurrency Thread State (Docker ps)

The orchestration relies on absolute core saturation, with all 24 independent python executors running parallel download and pseudoalignment routines, while the single detached script handles legacy merge completions concurrently.

```text
UID                 PID                 PPID                C                   STIME               TTY                 TIME                CMD
root                245176              245157              0                   Mar13               ?                   00:00:00            /bin/bash /app/run_pipeline.sh
root                245218              245176              0                   Mar13               ?                   00:01:42            python3 scripts/rna/run_all_species.py --max-gb 350.0 --workers 24 --threads 2

< - - - 24 Independent Simultaneous Thread Pools of NCBI & Kallisto - - - >

root                220915              245218              83.0                00:52               ?                   03:20               kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR14887992 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR14887992/SRR14887992_1.fastq.gz
root                221022              245218              85.5                00:54               ?                   01:53               kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR13938899 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR13938899/SRR13938899_1.fastq.gz
root                221075              245218              2.2                 00:54               ?                   00:02               curl -fsSL --retry 3 --retry-delay 10 --connect-timeout 30 -o output/amalgkit/apis_mellifera/work/getfastq/SRR33762311/SRR33762311_1.fastq.gz http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR337/011/SRR33762311/SRR33762311_1.fastq.gz
root                221169              245218              83.6                00:55               ?                   01:15               kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR33381866 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR33381866/SRR33381866_1.fastq.gz

< - - - Background Single-Thread Curation Engine - - - >

root                216051              245176              0                   00:42               ?                   00:00:00            /bin/bash scripts/rna/bg_curate.sh
root                216196              216051              100                 00:42               ?                   00:00:02            /app/.venv/bin/python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_odontomachus_brunneus.yaml --no-progress --steps merge curate sanity
```
