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
  apis_mellifera                 6284        4        0       20     1038       24    7370  ❌ Not run
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
  TOTAL                          6284        4        0       20     5106       38   11452
================================================================================
```

## 2. Progress Since Last Check

Since the previous check exactly 120 minutes ago (20:55 UTC to 22:55 UTC), the ThreadPool has smoothly quantified an additional **+235 `apis_mellifera` samples**, bringing the total completed for the honeybee from 803 to 1,038.

- **Running 10-Minute Average:** After the massive velocity spike two hours ago, the pipeline has comfortably settled into crunching larger SRA reads for `apis_mellifera` and is currently maintaining a highly consistent **19.6 samples quantified every 10 minutes** sustained.
- The threadpool has effectively consumed over 14% of the massive 7,370 Apis datasets queue.
- This continues to validate the pipeline's immense throughput bounds. The worker count remains saturated with 24 dedicated cores consistently grabbing work.

## 3. Live Thread Current Activities (Docker Top Snapshot)

The ThreadPool is maintaining 100% saturation and non-blocking concurrency:

```text
UID                 PID                 PPID                C                   STIME               TTY                 TIME                CMD
root                245176              245157              0                   Mar13               ?                   00:00:00            /bin/bash /app/run_pipeline.sh
root                245218              245176              0                   Mar13               ?                   00:01:41            python3 scripts/rna/run_all_species.py --max-gb 350.0 --workers 24 --threads 2
root                480194              245176              84                  20:26               ?                   02:05:07            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR29535102 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29535102/SRR29535102_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29535102/SRR29535102_2.fastq.gz
root                490437              245218              0                   21:52               ?                   00:00:02            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7283 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                490498              490437              78                  21:52               ?                   00:49:49            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR6031641 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR6031641/SRR6031641_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR6031641/SRR6031641_2.fastq.gz
root                494825              245218              0                   22:20               ?                   00:00:02            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7171 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                494869              494825              63                  22:20               ?                   00:22:17            kallisto quant --threads 1 --index output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR26053392 --single -l 200 -s 20.0 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR26053392/SRR26053392.fastq.gz
root                495198              245218              0                   22:22               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7181 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                495247              495198              72                  22:22               ?                   00:23:49            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR25389644 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR25389644/SRR25389644_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR25389644/SRR25389644_2.fastq.gz
root                495822              245218              0                   22:27               ?                   00:00:02            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 6626 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                495866              495822              73                  22:27               ?                   00:20:16            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR13871650 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR13871650/SRR13871650_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR13871650/SRR13871650_2.fastq.gz
root                496769              245218              0                   22:37               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7363 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                496812              496769              66                  22:37               ?                   00:12:17            kallisto quant --threads 1 --index output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR28818366 --single -l 200 -s 20.0 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR28818366/SRR28818366.fastq.gz
root                497143              245218              0                   22:40               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 6873 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                497190              497143              64                  22:40               ?                   00:09:47            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR15840297 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR15840297/SRR15840297_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR15840297/SRR15840297_2.fastq.gz
root                497581              245218              0                   22:43               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 4703 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                497623              497581              72                  22:43               ?                   00:08:43            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR37045114 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR37045114/SRR37045114_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR37045114/SRR37045114_2.fastq.gz
root                499074              245218              3                   22:54               ?                   00:00:02            curl -fsSL --retry 3 --retry-delay 10 --retry-connrefused --retry-all-errors --connect-timeout 30 -o output/amalgkit/apis_mellifera/work/getfastq/SRR26559880/SRR26559880_1.fastq.gz http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/080/SRR26559880/SRR26559880_1.fastq.gz
root                499119              499075              74                  22:54               ?                   00:00:59            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR29831676 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29831676/SRR29831676_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29831676/SRR29831676_2.fastq.gz
root                499190              499130              81                  22:54               ?                   00:00:35            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR20331285 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR20331285/SRR20331285_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR20331285/SRR20331285_2.fastq.gz
```
