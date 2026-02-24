import pandas as pd
df = pd.read_csv('output/amalgkit/acromyrmex_echinatior/work/metadata/metadata.tsv', sep='\t')
sizes = {}
for _, row in df.iterrows():
    try:
        spots = float(row.get('total_spots', 0))
        length = float(row.get('spot_length', 150))
        sizes[row['run']] = spots * length
    except Exception as e:
        print("err", e)
        sizes[row['run']] = 0.0
print("sum:", sum(sizes.values()))
