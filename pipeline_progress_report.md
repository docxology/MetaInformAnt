# Real-time Amalgkit Pipeline Progress

> [!TIP]
> **How to Monitor Pipeline Progress Yourself**
> You can pull this exact live telemetry natively from your local terminal at any time while the background orchestrator is running.
> Just run this command to safely check the SQL matrix and running threads via Docker:
> ```bash
> gcloud compute ssh metainformant-pipeline --zone=us-central1-a --project=cryptoptera --tunnel-through-iap --ssh-flag="-q" --command='sudo docker exec metainformant-pipeline-fresh python3 scripts/rna/check_pipeline_status.py -v && echo -e "\n--- PROCESS LIST ---\n" && sudo docker top metainformant-pipeline-fresh'
> ```

## 1. Live Progress Matrix

```text
================================================================================
  Amalgkit Pipeline Status  │  Process: 🟢 RUNNING
================================================================================
Species                       pending  downloa  downloa  quantif  quantif   failed   Total  Downstream      
--------------------------------------------------------------------------------
  apis_mellifera                 6166        9        0       15     1153       27    7370  ❌ Not run
  harpegnathos_saltator             0        0        0        0      689        0     689  ✅ Complete
  temnothorax_longispinosus         0        0        0        0      508        0     508  ⚠️  Merge only
  solenopsis_invicta                0        0        0        0      450        1     451  ⚠️  Merge only
  monomorium_pharaonis              0        0        0        0      370        0     370  ⚠️  Merge only
  camponotus_floridanus             0        0        0        0      366        1     367  ⚠️  Merge only
  temnothorax_americanus            0        0        0        0      331        0     331  ❌ Not run
  nasonia_vitripennis               0        0        0        0      280        0     280  ❌ Not run
  ooceraea_biroi                    0        0        0        0      274        0     274  ❌ Not run
  helicoverpa_armigera              0        0        0        0      239        0     239  ❌ Not run
  atta_cephalotes                   0        0        0        0      217        0     217  ✅ Complete
  athalia_rosae                     0        0        0        0      197        0     197  ❌ Not run
  bombus_terrestris                 0        0        0        0      180        0     180  ❌ Not run
  linepithema_humile                0        0        0        0      173        0     173  ❌ Not run
  cardiocondyla_obscurior           0        0        0        0      167        0     167  ⚠️  Merge only
  temnothorax_nylanderi             0        0        0        0      154       12     166  ⚠️  Merge only
  spodoptera_frugiperda             0        0        0        0      143        0     143  ❌ Not run
  pogonomyrmex_barbatus             0        0        0        0      132        0     132  ❌ Not run
  blattella_germanica               0        0        0        0      101        0     101  ❌ Not run
  megachile_rotundata               0        0        0        0       95        0      95  ❌ Not run
  polistes_fuscatus                 0        0        0        0       60        0      60  ❌ Not run
  acromyrmex_echinatior             0        0        0        0       44        0      44  ✅ Complete
  temnothorax_curvispinosus         0        0        0        0       43        0      43  ⚠️  Merge only
  polistes_canadensis               0        0        0        0       41        0      41  ❌ Not run
  nylanderia_fulva                  0        0        0        0       40        0      40  ⚠️  Merge only
  trichogramma_pretiosum            0        0        0        0       40        0      40  ❌ Not run
  wasmannia_auropunctata            0        0        0        0       33        0      33  ❌ Not run
  coptotermes_formosanus            0        0        0        0       29        0      29  ❌ Not run
  formica_exsecta                   0        0        0        0       23        0      23  ⚠️  Merge only
  odontomachus_brunneus             0        0        0        0       19        0      19  ⚠️  Merge only
  vollenhovia_emeryi                0        0        0        0       15        0      15  ❌ Not run
  dinoponera_quadriceps             0        0        0        0       13        0      13  ❌ Not run
  anoplolepis_gracilipes            0        0        0        0        7        0       7  ❌ Not run
--------------------------------------------------------------------------------
  TOTAL                          6166        9        0       15     5221       41   11452
================================================================================
```

## 2. Progress Since Last Check

Since the previous check exactly 32 minutes ago (23:13 UTC to 23:45 UTC), the ThreadPool has successfully quantified an additional **+70 `apis_mellifera` samples**, bringing the total completed for the honeybee from 1,083 to 1,153.

- **Running 10-Minute Average:** The pipeline velocity remains exceptionally steady at **21.8 samples quantified every 10 minutes** over this half-hour window.
- The 24 worker threads remain exclusively locked onto the 7,370 Apis datsets queue.
- **Post-Merge Curations:** The 10 distinct species marked `⚠️ Merge only` have absolutely **no background processes** actively running downstream components (e.g. `cstmm`, `curate`). All 24 threads are exclusively dumping and quantifying `apis_mellifera` FASTQ files.

## 3. Live Thread Current Activities (Docker Top Snapshot)

The ThreadPool is maintaining 100% saturation and non-blocking concurrency:

```text
UID                 PID                 PPID                C                   STIME               TTY                 TIME                CMD
root                245176              245157              0                   Mar13               ?                   00:00:00            /bin/bash /app/run_pipeline.sh
root                245218              245176              0                   Mar13               ?                   00:01:42            python3 scripts/rna/run_all_species.py --max-gb 350.0 --workers 24 --threads 2
root                480194              245176              82                  20:26               ?                   02:44:36            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR29535102 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29535102/SRR29535102_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29535102/SRR29535102_2.fastq.gz
root                496769              245218              0                   22:37               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7363 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                496812              496769              70                  22:37               ?                   00:49:10            kallisto quant --threads 1 --index output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR28818366 --single -l 200 -s 20.0 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR28818366/SRR28818366.fastq.gz
root                502395              245218              0                   23:16               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7068 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                502464              502395              75                  23:17               ?                   00:22:42            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR15560394 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR15560394/SRR15560394_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR15560394/SRR15560394_2.fastq.gz
root                503865              245218              0                   23:26               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7076 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                503909              503865              67                  23:27               ?                   00:13:35            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR25389648 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR25389648/SRR25389648_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR25389648/SRR25389648_2.fastq.gz
root                504010              245218              0                   23:27               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 6685 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                504053              504010              80                  23:27               ?                   00:15:31            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR13871636 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR13871636/SRR13871636_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR13871636/SRR13871636_2.fastq.gz
root                504670              245218              0                   23:33               ?                   00:00:02            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 6707 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                504712              504670              74                  23:33               ?                   00:10:17            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR23183309 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR23183309/SRR23183309_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR23183309/SRR23183309_2.fastq.gz
root                504713              245218              0                   23:33               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7238 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                504737              245218              0                   23:33               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 6759 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                504796              504713              62                  23:33               ?                   00:08:30            kallisto quant --threads 1 --index output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR2954346 --single -l 200 -s 20.0 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR2954346/SRR2954346.fastq.gz
root                504798              504737              73                  23:33               ?                   00:09:59            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR23183246 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR23183246/SRR23183246_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR23183246/SRR23183246_2.fastq.gz
root                504883              245218              0                   23:34               ?                   00:00:02            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 6515 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                504925              504883              69                  23:34               ?                   00:08:46            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR15840327 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR15840327/SRR15840327_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR15840327/SRR15840327_2.fastq.gz
root                505367              245218              99                  23:38               ?                   00:10:47            pigz -p 2 output/amalgkit/apis_mellifera/work/getfastq/DRR493329/DRR493329_1.fastq
```
