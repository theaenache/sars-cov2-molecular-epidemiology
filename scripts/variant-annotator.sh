#!/bin/bash

# variant-annotator.sh
# variant annotation time
# annotates SARS-CoV-2 variants using Ensembl VEP

# 1) comment out header lines in covid19-snps.tsv
sed -i 's/^\//#\//' covid19-snps.tsv
sed -i 's/^NUCMER/#NUCMER/' covid19-snps.tsv
sed -i 's/^\[P1\]/#[P1]/' covid19-snps.tsv
sed -i '/^$/d' covid19-snps.tsv

# 2) convert MUMmer SNP file to VCF format
python mummer-snps2vcf covid19-snps.tsv --reference reference-covid-19.fna > covid19-snps.vcf

# 3) add contig header line required by VEP
sed -i '2i ##contig=<ID=NC_045512.2,length=29903>' covid19-snps.vcf

# 4) create and activate conda environment
conda create -n vep-env -y
source /opt/conda/etc/profile.d/conda.sh
conda activate vep-env

# 5) install VEP and dependencies
mamba install bioconda::ensembl-vep biopython -y

# 6) index the GTF annotation file
tabix -p gff reference-covid-19-4vep.gtf.gz

# 7) run VEP annotation
nohup vep -i covid19-snps.vcf \
    --gtf reference-covid-19-4vep.gtf.gz \
    --fasta reference-covid-19.fna \
    --vcf \
    -o covid19-snps-annotated.vcf \
    --force_overwrite > vep.log 2>&1

# 8) filter VEP output by impact
python vep-filter-impact covid19-snps-annotated.vcf > covid19-snps-annotated-filtered.vcf

echo "the annotated variants were written to covid19-snps-annotated-filtered.vcf"
