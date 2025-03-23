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

3. Merge AFR, AMR, EUR, FIN data
```bash
bash merge.bsh
```

4. Principal Component Analysis (PCA) </br>
```bash
bash pca.bsh

# based on the scree plots we decided to use 5 PCs for each ancestry for association analysis
Rscript scree_plot.R

# plot PCs by affection status and sex for each ancestry
Rscript mega_pca_aff_${anc}.R
Rscript mega_pca_sex_${anc}.R
```
