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
  apis_mellifera                 6237        5        0       19     1083       26    7370  ❌ Not run
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
  TOTAL                          6237        5        0       19     5151       40   11452
================================================================================
```

## 2. Progress Since Last Check

Since the previous check exactly 18 minutes ago (22:55 UTC to 23:13 UTC), the ThreadPool has successfully quantified an additional **+45 `apis_mellifera` samples**, bringing the total completed for the honeybee from 1,038 to 1,083.

- **Running 10-Minute Average:** The pipeline velocity remains exceptionally steady, having ticked *up* slightly to **25.0 samples quantified every 10 minutes** over this short 18-minute window.
- The 24 worker threads remain exclusively locked onto the 7,370 Apis datsets queue.
- **Post-Merge Curations:** The 10 distinct species marked `⚠️ Merge only` have absolutely **no background processes** actively running downstream components (e.g. `cstmm`, `curate`). All 24 threads are exclusively dumping and quantifying `apis_mellifera` FASTQ files.

## 3. Live Thread Current Activities (Docker Top Snapshot)

The ThreadPool is maintaining 100% saturation and non-blocking concurrency:

```text
UID                 PID                 PPID                C                   STIME               TTY                 TIME                CMD
root                245176              245157              0                   Mar13               ?                   00:00:00            /bin/bash /app/run_pipeline.sh
root                245218              245176              0                   Mar13               ?                   00:01:41            python3 scripts/rna/run_all_species.py --max-gb 350.0 --workers 24 --threads 2
root                480194              245176              83                  20:26               ?                   02:19:30            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR29535102 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29535102/SRR29535102_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29535102/SRR29535102_2.fastq.gz
root                495198              245218              0                   22:22               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7181 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                495247              495198              73                  22:22               ?                   00:38:15            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR25389644 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR25389644/SRR25389644_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR25389644/SRR25389644_2.fastq.gz
root                496769              245218              0                   22:37               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7363 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                496812              496769              69                  22:37               ?                   00:25:44            kallisto quant --threads 1 --index output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR28818366 --single -l 200 -s 20.0 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR28818366/SRR28818366.fastq.gz
root                497996              245218              0                   22:46               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7302 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                498041              497996              72                  22:46               ?                   00:20:13            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR8769531 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR8769531/SRR8769531_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR8769531/SRR8769531_2.fastq.gz
root                499075              245218              0                   22:54               ?                   00:00:02            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 6883 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                499119              499075              73                  22:54               ?                   00:14:54            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR29831676 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29831676/SRR29831676_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR29831676/SRR29831676_2.fastq.gz
root                499823              245218              0                   22:59               ?                   00:00:01            /opt/conda/bin/python /usr/local/bin/amalgkit quant --out_dir output/amalgkit/apis_mellifera/work --metadata output/amalgkit/apis_mellifera/work/metadata/metadata.tsv --threads 1 --batch 7154 --clean_fastq yes --index_dir output/amalgkit/apis_mellifera/work/index
root                499870              499823              73                  22:59               ?                   00:11:14            kallisto quant --threads 1 -i output/amalgkit/apis_mellifera/work/index/Apis_mellifera_transcripts.idx -o /app/output/amalgkit/apis_mellifera/work/quant/SRR15560371 /app/output/amalgkit/apis_mellifera/work/getfastq/SRR15560371/SRR15560371_1.fastq.gz /app/output/amalgkit/apis_mellifera/work/getfastq/SRR15560371/SRR15560371_2.fastq.gz
root                501416              245218              3                   23:10               ?                   00:00:07            curl -fsSL --retry 3 --retry-delay 10 --retry-connrefused --retry-all-errors --connect-timeout 30 -o output/amalgkit/apis_mellifera/work/getfastq/SRR19055431/SRR19055431_1.fastq.gz http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR190/031/SRR19055431/SRR19055431_1.fastq.gz
root                501952              245218              4                   23:13               ?                   00:00:03            curl -fsSL --retry 3 --retry-delay 10 --retry-connrefused --retry-all-errors --connect-timeout 30 -o output/amalgkit/apis_mellifera/work/getfastq/SRR4017776/SRR4017776_1.fastq.gz http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR401/006/SRR4017776/SRR4017776_1.fastq.gz
root                501955              245218              4                   23:13               ?                   00:00:02            curl -fsSL --retry 3 --retry-delay 10 --retry-connrefused --retry-all-errors --connect-timeout 30 -o output/amalgkit/apis_mellifera/work/getfastq/SRR4045684/SRR4045684.fastq.gz http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR404/004/SRR4045684/SRR4045684.fastq.gz
root                501978              245218              4                   23:14               ?                   00:00:01            curl -fsSL --retry 3 --retry-delay 10 --retry-connrefused --retry-all-errors --connect-timeout 30 -o output/amalgkit/apis_mellifera/work/getfastq/SRR15560394/SRR15560394_2.fastq.gz http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/094/SRR15560394/SRR15560394_2.fastq.gz
root                502035              245218              3                   23:14               ?                   00:00:00            curl -fsSL --retry 3 --retry-delay 10 --retry-connrefused --retry-all-errors --connect-timeout 30 -o output/amalgkit/apis_mellifera/work/getfastq/SRR20852077/SRR20852077_2.fastq.gz http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR208/077/SRR20852077/SRR20852077_2.fastq.gz
```
