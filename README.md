# SARS-CoV-2 Molecular Epidemiology Pipeline
*Computational molecular epidemiology pipeline for SARS-CoV-2 variant analysis*



---

## Overview

This project investigates genetic variants associated with virulence in SARS-CoV-2 using whole-genome sequencing data from 272,726 clinical isolates collected in the United Kingdom. The pipeline covers variant calling, genome-wide association analysis, variant annotation, and phylogenetic tree construction.

We found that variant **3653_C_T** in **ORF1ab** (L1130F — Leucine → Phenylalanine at position 1130) was significantly associated with higher viral load (lower cycle threshold), with a mean CT of 21.34 among carriers vs. 34.03 among non-carriers (p < 2.2×10⁻¹⁶).

---

## Data

| File | Description |
|------|-------------|
| `covid19-genomes-PRJEB37886-dec2021.fna.gz` | 272,726 assembled SARS-CoV-2 genomes (NCBI BioProject PRJEB37886) |
| `reference-covid-19.fna` | Wuhan-Hu-1 reference genome (NC_045512) |
| `reference-covid-19-4vep.gtf.gz` | Reference genome annotation (VEP-compatible) |
| `reference-covid-19.gb` | Reference genome in GenBank format |
| `covid19_cycle_threshold_mock_data.csv` | Mock PCR cycle threshold values per sample |

Samples were sequenced on the Illumina NovaSeq 6000 platform by the COVID-19 Genomics UK (COG-UK) consortium between November 29, 2021 and January 15, 2022.

---

## Pipeline

### Part 1 — Variant Calling
**Script:** `scripts/call-variants.sh`

Aligns all 272,726 SARS-CoV-2 genomes against the Wuhan-Hu-1 reference using MUMmer's `dnadiff` wrapper and extracts SNPs with `show-snps`.

```bash
bash scripts/call-variants.sh
```

**Output:** `covid19-snps.tsv.gz` (992MB uncompressed, 26,726 unique variants)

---

### Part 2 — Genome-Wide Association with Virulence
**Script:** `scripts/gwa-ct.py`

Tests each variant for association with PCR cycle threshold (CT) values using a Student's t-test. Lower CT = higher viral load = higher virulence.

```bash
python3 scripts/gwa-ct.py
```

**Output:** `results/covid19_vars_ct.csv`

| Column | Description |
|--------|-------------|
| Variant | Position and nucleotide change (e.g. 3653_C_T) |
| Samples | Number of samples carrying the variant |
| Mean_CT | Mean cycle threshold among carriers |
| SD_CT | Standard deviation of CT among carriers |
| t_test | p-value from Student's t-test |

---

### Part 3 — Variant Annotation
**Script:** `scripts/variant-annotator.sh`

Annotates variants using Ensembl VEP with a SARS-CoV-2-specific GTF annotation file.

```bash
bash scripts/variant-annotator.sh
```

**Output:** `results/covid19-snps-annotated-filtered.vcf.gz`

---

### Part 4 — Phylogenetic Tree
**Scripts:** `scripts/snp-tree.sh`, `scripts/build_phylip.py`

Builds a maximum likelihood phylogenetic tree from 500 representative SARS-CoV-2 samples using SARS-CoV-1 (NC_004718) as an outgroup. Best substitution model (TPM1uf+G4) selected by modeltest-ng.

```bash
bash scripts/snp-tree.sh
```

**Output:** `results/covid_tree.raxml.bestTree` (visualized on iTOL)

---

## Results Summary

| Metric | Value |
|--------|-------|
| Total samples | 272,726 |
| Unique variants | 26,726 |
| Top variant | 3653_C_T |
| Gene | ORF1ab |
| Amino acid change | L1130F (Leucine → Phenylalanine) |
| Samples with variant | 1,357 |
| Mean CT (with variant) | 21.34 |
| Mean CT (without variant) | 34.03 |
| p-value | < 2.2×10⁻¹⁶ |

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| MUMmer | 3.23 | Whole genome alignment |
| Python | 3.12 | GWAS, phylip matrix construction |
| scipy | — | Student's t-test |
| Ensembl VEP | — | Variant annotation |
| bcftools | 1.23 | VCF manipulation |
| modeltest-ng | 0.1.7 | Substitution model selection |
| raxml-ng | — | Maximum likelihood tree |

---

## Repository Structure

```
├── scripts/
│   ├── call-variants.sh       # Part 1: variant calling
│   ├── gwa-ct.py              # Part 2: GWAS
│   ├── variant-annotator.sh   # Part 3: VEP annotation
│   ├── snp-tree.sh            # Part 4: phylogenetic tree
│   └── build_phylip.py        # builds PHYLIP matrix from TSV
├── results/
│   ├── covid19_vars_ct.csv               # GWAS results
│   ├── covid19-snps.tsv.gz               # variant calls (compressed)
│   ├── covid19-snps-annotated-filtered.vcf.gz  # annotated variants
│   └── covid_tree.raxml.bestTree         # best ML tree
├── writeup/
│   └── writeup.pdf            # project writeup
└── README.md
```
