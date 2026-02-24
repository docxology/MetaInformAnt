import os
import glob
import pandas as pd
import time

print("| Species | Total | Downloaded | Quantified | Rate | ETA (Samples) | ETA (Data) |")
print("|---|---|---|---|---|---|---|")

start_time = None
pid_file = 'output/amalgkit/run_all_species_incremental.pid'
if os.path.exists(pid_file):
    start_time = os.path.getmtime(pid_file)

config_files = sorted(glob.glob('config/amalgkit/amalgkit_*.yaml'))
for cfg in config_files:
    if 'template' in cfg or '_test.yaml' in cfg or 'cross_species' in cfg:
        continue
    
    species = os.path.basename(cfg).replace('amalgkit_', '').replace('.yaml', '')
    metadata_file = f'output/amalgkit/{species}/work/metadata/metadata.tsv'
    
    try:
        df = pd.read_csv(metadata_file, sep='\t')
        total_samples = len(df)
        runs = df['run'].tolist()
        
        # approximate data sizes from spots and length
        sizes = {}
        for _, row in df.iterrows():
            try:
                spots = float(row.get('total_spots', 0))
                length = float(row.get('spot_length', 150))
                if pd.isna(spots): spots = 0.0
                if pd.isna(length): length = 150.0
                sizes[row['run']] = spots * length
            except:
                sizes[row['run']] = 0.0
    except Exception as e:
        total_samples = 0
        runs = []
        sizes = {}

    downloaded = 0
    quantified = 0
    quantified_size = 0.0
    total_size = sum(sizes.values())

    getfastq_dir = f'output/amalgkit/{species}/fastq/getfastq'
    quant_dir = f'output/amalgkit/{species}/work/quant'

    for run in runs:
        is_quantified = False
        if os.path.exists(f'{quant_dir}/{run}/{run}_abundance.tsv'):
            is_quantified = True
            quantified += 1
            quantified_size += sizes.get(run, 0.0)
            
        is_downloaded = False
        run_getfastq = f'{getfastq_dir}/{run}'
        if is_quantified is True:
            is_downloaded = True
        elif os.path.isdir(run_getfastq) and os.listdir(run_getfastq):
            is_downloaded = True
            
        if is_downloaded:
            downloaded += 1

    eta_samples = "N/A"
    eta_data = "N/A"
    rate_str = "N/A"
    
    if start_time and quantified > 0 and total_samples > 0:
        time_elapsed = time.time() - start_time
        if time_elapsed > 0:
            samples_rate = quantified / time_elapsed
            rate_str = f"{samples_rate * 3600:.1f} smp/h"
            
            remaining_samples = total_samples - quantified
            if remaining_samples > 0:
                eta_sec = remaining_samples / samples_rate
                eta_samples = f"{eta_sec/3600:.1f}h"
            else:
                eta_samples = "Done"
                
            if total_size > 0:
                data_rate = quantified_size / time_elapsed
                remaining_size = total_size - quantified_size
                if data_rate > 0 and remaining_size > 0:
                    eta_sec = remaining_size / data_rate
                    eta_data = f"{eta_sec/3600:.1f}h"
                elif remaining_size <= 0:
                    eta_data = "Done"

    print(f"| {species} | {total_samples} | {downloaded} | {quantified} | {rate_str} | {eta_samples} | {eta_data} |")
