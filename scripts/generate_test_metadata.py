import csv
from pathlib import Path

path = Path("output/amalgkit/pbarbatus_test5/work/metadata/metadata_selected.tsv")
path.parent.mkdir(parents=True, exist_ok=True)

data = [
    {"scientific_name": "Pogonomyrmex barbatus", "sample_group": "ovary", "run": "SRR34065669", "AWS_Link": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR34065669/SRR34065669"},
    {"scientific_name": "Pogonomyrmex barbatus", "sample_group": "brain", "run": "SRR14740521", "AWS_Link": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14740521/SRR14740521"},
    {"scientific_name": "Pogonomyrmex barbatus", "sample_group": "brain", "run": "SRR1817181", "AWS_Link": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1817181/SRR1817181"},
    {"scientific_name": "Pogonomyrmex barbatus", "sample_group": "brain", "run": "SRR1817181_DUP", "AWS_Link": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1817181/SRR1817181"},
    {"scientific_name": "Pogonomyrmex barbatus", "sample_group": "ovary", "run": "SRR34065669_DUP", "AWS_Link": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR34065669/SRR34065669"},
]

with open(path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["scientific_name", "sample_group", "run", "AWS_Link"], delimiter="\t")
    writer.writeheader()
    writer.writerows(data)

print("Metadata file created successfully.")
