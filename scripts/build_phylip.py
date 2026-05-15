import csv
from collections import defaultdict

with open("samples500.txt") as f:
    samples = [line.strip() for line in f if line.strip()]
sample_set = set(samples)

print("Parsing covid SNPs...")
variant_samples = defaultdict(set)
with open("covid19-snps.tsv") as f:
    for line in f:
        if line.startswith("#") or line.strip() == "":
            continue
        cols = line.strip().split("\t")
        if len(cols) < 12:
            continue
        pos, ref, alt, sample = cols[0], cols[1], cols[2], cols[11]
        if ref == "." or alt == ".":
            continue
        if sample in sample_set:
            variant_samples[f"{pos}_{ref}_{alt}"].add(sample)

print("Parsing outgroup SNPs...")
outgroup_variants = set()
with open("outgroup-snps.tsv") as f:
    for line in f:
        if line.startswith("#") or line.strip() == "":
            continue
        cols = line.strip().split("\t")
        if len(cols) < 12:
            continue
        pos, ref, alt = cols[0], cols[1], cols[2]
        if ref == "." or alt == ".":
            continue
        outgroup_variants.add(f"{pos}_{ref}_{alt}")

all_variants = sorted(variant_samples.keys() | outgroup_variants,
                      key=lambda x: int(x.split("_")[0]))
all_samples = samples + ["NC_004718_outgroup"]

print(f"Unique variants: {len(all_variants)}")
print(f"Total samples: {len(all_samples)}")

with open("covid_tree.phy", "w") as out:
    out.write(f" {len(all_samples)} {len(all_variants)}\n")
    for sample in samples:
        seq = "".join("1" if sample in variant_samples[v] else "0"
                      for v in all_variants)
        out.write(f"{sample[:10]:<10} {seq}\n")
    seq = "".join("1" if v in outgroup_variants else "0"
                  for v in all_variants)
    out.write(f"{'outgroup':<10} {seq}\n")

print("Done! Written to covid_tree.phy")
