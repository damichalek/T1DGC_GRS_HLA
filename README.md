# T1DGC: HLA-focused type 1 diabetes genetic risk prediction in populations of diverse ancestry

This repository contains the steps taken for HLA-focused type 1 diabetes genetic risk prediction in Admixed African (AFR), Admixed American (AMR), European (EUR) and Finnish (FIN) ancestry individuals. </br>

For questions, please contact Dominika Michalek (dam8mt@virginia.edu).

## Genotyped data
DNA samples were genotyped on the Illumina ImmunoChip array. </br>

Genotyped data underwent following steps:
1. Quality control </br>

Genotype data in the HLA region on human chromosome 6 (28Mbp - 34Mbp) were extracted from the ImmunoChip panel, with quality control performed as previously published ([Robertson CC, Inshaw JRJ, et al 2020](https://pubmed.ncbi.nlm.nih.gov/34127860/)).

2. Population structure and ancestry inference
```bash
bash anc_infer.bsh
```
