import csv
from scipy import stats
import numpy as np

# load our cycle threshold data
ct_data = {}
with open("covid19_cycle_threshold_mock_data.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        ct_data[row["Isolate"].strip()] = float(row["Cycle_Threshold"])

# parse through covid19-snps.tsv and collect samples per variant
variant_samples = {}
with open("covid19-snps.tsv") as f:
    for line in f:
        if line.startswith("#") or line.startswith("/") or line.startswith("NUCMER") or line.strip() == "" or line.startswith("["):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 12:
            continue
        pos = cols[0]
        ref_base = cols[1]
        alt_base = cols[2]
        sample = cols[11]
        if ref_base == "." or alt_base == ".":
            continue
        variant = f"{pos}_{ref_base}_{alt_base}"
        if variant not in variant_samples:
            variant_samples[variant] = set()
        variant_samples[variant].add(sample)

# test each variant for association with CT
results = []
for variant, samples in variant_samples.items():
    ct_with = [ct_data[s] for s in samples if s in ct_data]
    if len(ct_with) < 2:
        continue
    ct_without = [ct_data[s] for s in ct_data if s not in samples]
    mean_ct = np.mean(ct_with)
    sd_ct = np.std(ct_with, ddof=1)
    t_stat, p_val = stats.ttest_ind(ct_with, ct_without)
    results.append({
        "Variant": variant,
        "Samples": len(ct_with),
        "Mean_CT": round(mean_ct, 4),
        "SD_CT": round(sd_ct, 4),
        "t_test": p_val
    })

# write our output
with open("covid19_vars_ct.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Variant", "Samples", "Mean_CT", "SD_CT", "t_test"])
    writer.writeheader()
    writer.writerows(results)

print(f"done!! tested {len(results)} variants. output written to covid19_vars_ct.csv")