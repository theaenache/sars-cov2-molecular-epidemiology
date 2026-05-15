#!/bin/bash

# snp-tree.sh
# the phylogenetic tree time
# builds a phylogenetic tree from SARS-CoV-2 SNP data
# using a SARS-CoV-1 outgroup (NC_004718)


# 1) download SARS-CoV-1 outgroup sequence from NCBI
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_004718&rettype=fasta&retmode=text" -O outgroup.fna

# 2) call variants between outgroup and reference using MUMmer
dnadiff reference-covid-19.fna outgroup.fna -p outgroup
show-snps -rT outgroup.delta > outgroup-snps.tsv

# 3) comment out header lines in outgroup SNP file
sed -i 's/^\//#\//' outgroup-snps.tsv
sed -i 's/^NUCMER/#NUCMER/' outgroup-snps.tsv
sed -i 's/^\[P1\]/#[P1]/' outgroup-snps.tsv
sed -i '/^$/d' outgroup-snps.tsv

# 4) select 500 random samples from covid SNP file
grep -v "^#" covid19-snps.tsv | cut -f12 | sort -u | shuf -n 500 > samples500.txt

# 5) build SNP presence/absence matrix in PHYLIP format
# using a custom Python script (build_phylip.py)
python3 build_phylip.py

# 6) convert binary characters to nucleotides for modeltest-ng
python3 -c "
with open('covid_tree.phy') as f:
    lines = f.readlines()
with open('covid_tree.phy', 'w') as out:
    out.write(lines[0])
    for line in lines[1:]:
        name = line[:11]
        seq = line[11:].strip()
        seq = seq.replace('1', 'T').replace('0', 'A')
        out.write(f'{name}{seq}\n')
"

# 7) install modeltest-ng and raxml-ng
mamba install bioconda::raxml-ng bioconda::modeltest-ng -y

# 8) run modeltest-ng to find best substitution model
modeltest-ng -i covid_tree.phy -d nt -o modeltest_out

# 9) run raxml-ng to build phylogenetic tree
# out best model is determined by modeltest-ng: TPM1uf+G4
raxml-ng --msa covid_tree.phy \
    --model TPM1uf+G4 \
    --prefix covid_tree \
    --outgroup outgroup \
    --seed 12345

# 10) best tree is in covid_tree.raxml.bestTree
# then, we upload to iTOL at https://itol.embl.de/upload.cgi for visualization
echo "best tree written to covid_tree.raxml.bestTree"
echo "now, upload covid_tree.raxml.bestTree to https://itol.embl.de/upload.cgi !"
