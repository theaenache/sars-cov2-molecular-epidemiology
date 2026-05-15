#!/bin/bash

gunzip -c /home/jovyan/shared/ph-671/final_project_inputs/covid19-genomes-PRJEB37886-dec2021.fna.gz > covid19-genomes.fna

dnadiff reference-covid-19.fna covid19-genomes.fna

show-snps -rT out.delta > covid19-snps.tsv
