#This script must be run in the same directory as the score_only folder (from ProteinMPNN)
#To run in terminal: python3 custom_runs/summarize_scores.py
import pandas as pd
import re

score_file = "../outputs/my_variants/5UOI_sample_5_score_only_from_fasta/score_only/5UOI_fasta_all.txt"
data = []
with open(score_file, "r") as f:
    for line in f:
        m = re.match(
            r"Score for ([^ ]+) (?:from PDB|from FASTA), mean: ([\d\.]+), std: ([\d\.]+), sample size: (\d+),\s+global score, mean: ([\d\.]+), std: ([\d\.]+), sample size: (\d+)",
            line.strip()
        )
        if m:
            seq_name = m.group(1)
            mean_score = float(m.group(2))
            std_score = float(m.group(3))
            sample_size = int(m.group(4))
            global_mean = float(m.group(5))
            global_std = float(m.group(6))
            global_sample_size = int(m.group(7))
            data.append({
                "Sequence": seq_name,
                "Mean Score": mean_score,
                "Std Dev": std_score,
                "Sample Size": sample_size,
                "Global Mean": global_mean,
                "Global Std": global_std,
                "Global Sample Size": global_sample_size
            })
df = pd.DataFrame(data)
df = df.sort_values("Mean Score")
print(df.to_string(index=False))
df.to_csv("../outputs/my_variants/5UOI_sample_5_score_only_from_fasta/score_only/score_summary.csv", index=False)