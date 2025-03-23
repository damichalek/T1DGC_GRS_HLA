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
bash pca_merged.bsh

# based on the scree plots we decided to use 5 PCs for each ancestry and merged ancestry data for association analysis
Rscript scree_plot.R
Rscript scree_plot_merged.R

# plot PCs by affection status and sex for each ancestry and merged ancestry data
Rscript mega_pca_aff_${anc}.R
Rscript mega_pca_sex_${anc}.R
Rscript mega_pca_aff_merged.R
Rscript mega_pca_sex_merged.R
```

## HLA imputation
1. Impute to Multi-ethnic HLA reference panel (Minimac4) 1.7.3 on [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!).
- Reference Panel - four-digit multi-ethnic HLA reference panel v2 (GRCh37/hg19)
- Array Build - GRCh37/hg19
- Phasing - Eagle v2.4 (phased output)
- Mode - Quality Control & Imputation

2. Filter HLA imputation results (HLA alleles and HLA genes (MAF > 0.005 and R2 > 0.5))
```bash
bash HLA_filter_imp.bsh
bash HLA_filter_imp_AMR.bsh # after SAS separation from AMR
```

3. Remove individuals with missing phenotype information
```bash
bash HLA_remove_missing_pheno.bsh
```

4. Merge filtered imputed AFR, AMR, EUR, FIN data
```bash
bash HLA_merge.bsh
```

## HLA association analysis
