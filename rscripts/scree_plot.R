library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("M:/T1DGC/USERS/dam8mt/data/MEGA")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eigenval_AFR = read.delim("PCA/mega_pca_b37_pruned_unrelated_AFR_controls.eigenval", header = F)
eigenval_AMR = read.delim("PCA/mega_pca_b37_pruned_unrelated_AMR_controls.eigenval", header = F)
eigenval_AMR_old = read.delim("PCA/AMR_old/mega_pca_b37_pruned_unrelated_AMR_controls.eigenval", header = F)
eigenval_EAS = read.delim("PCA/mega_pca_b37_pruned_unrelated_EAS_controls.eigenval", header = F)
eigenval_EUR = read.delim("PCA/mega_pca_b37_pruned_unrelated_EUR_controls.eigenval", header = F)
eigenval_FIN = read.delim("PCA/mega_pca_b37_pruned_unrelated_FIN_controls.eigenval", header = F)
eigenval_SAS = read.delim("PCA/mega_pca_b37_pruned_unrelated_SAS_controls.eigenval", header = F)

# Convert column into rows
eigenval_AFR = as.data.frame(t(eigenval_AFR))
eigenval_AMR = as.data.frame(t(eigenval_AMR))
eigenval_AMR_old = as.data.frame(t(eigenval_AMR_old))
eigenval_EAS = as.data.frame(t(eigenval_EAS))
eigenval_EUR = as.data.frame(t(eigenval_EUR))
eigenval_FIN = as.data.frame(t(eigenval_FIN))
eigenval_SAS = as.data.frame(t(eigenval_SAS))

# Name the columns
colnames(eigenval_AFR) = c(paste0("PC", 1:10))
colnames(eigenval_AMR) = c(paste0("PC", 1:10))
colnames(eigenval_AMR_old) = c(paste0("PC", 1:10))
colnames(eigenval_EAS) = c(paste0("PC", 1:10))
colnames(eigenval_EUR) = c(paste0("PC", 1:10))
colnames(eigenval_FIN) = c(paste0("PC", 1:10))
colnames(eigenval_SAS) = c(paste0("PC", 1:10))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data frame
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Sum of eigenvalues (PC1:PC10)
sum_AFR = rowSums(eigenval_AFR)
sum_AMR = rowSums(eigenval_AMR)
sum_AMR_old = rowSums(eigenval_AMR_old)
sum_EAS = rowSums(eigenval_EAS)
sum_EUR = rowSums(eigenval_EUR)
sum_FIN = rowSums(eigenval_FIN)
sum_SAS = rowSums(eigenval_SAS)

# Proportion of variance - divide each eigenvalue (PC1:PC10) by sum of eigenvalues
var_explained_AFR = mapply('/', eigenval_AFR, sum_AFR)
var_explained_AMR = mapply('/', eigenval_AMR, sum_AMR)
var_explained_AMR_old = mapply('/', eigenval_AMR_old, sum_AMR_old)
var_explained_EAS = mapply('/', eigenval_EAS, sum_EAS)
var_explained_EUR = mapply('/', eigenval_EUR, sum_EUR)
var_explained_FIN = mapply('/', eigenval_FIN, sum_FIN)
var_explained_SAS = mapply('/', eigenval_SAS, sum_SAS)

# A data frame containing the PCs and the variance explained by each PC
var_explained_AFR_df = data.frame(PC=(1:10), var_explained=var_explained_AFR)
var_explained_AMR_df = data.frame(PC=(1:10), var_explained=var_explained_AMR)
var_explained_AMR_old_df = data.frame(PC=(1:10), var_explained=var_explained_AMR_old)
var_explained_EAS_df = data.frame(PC=(1:10), var_explained=var_explained_EAS)
var_explained_EUR_df = data.frame(PC=(1:10), var_explained=var_explained_EUR)
var_explained_FIN_df = data.frame(PC=(1:10), var_explained=var_explained_FIN)
var_explained_SAS_df = data.frame(PC=(1:10), var_explained=var_explained_SAS)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Scree plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

screePlot = function(data, title, path, filename){
  data %>% 
    ggplot(aes(x = PC, y = var_explained, group = 1)) +
    scale_x_continuous(n.breaks = 10) +
    geom_point(size = 2) +
    geom_line() +
    ylab("Variance explained") +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0(path,filename), width=15, height=15, units="cm", dpi=300)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dataList = list(var_explained_AFR_df, var_explained_AMR_df, var_explained_AMR_old_df, var_explained_EAS_df, var_explained_EUR_df, var_explained_FIN_df, var_explained_SAS_df)
titleList = list("Scree plot AFR", "Scree plot AMR", "Scree plot AMR old", "Scree plot EAS", "Scree plot EUR", "Scree plot FIN", "Scree plot SAS")
filenameList = list("scree_plot_AFR.png", "scree_plot_AMR.png", "scree_plot_AMR_old.png", "scree_plot_EAS.png", "scree_plot_EUR.png","scree_plot_FIN.png", "scree_plot_SAS.png")

mapply(screePlot, dataList, titleList, "Graphs/PCA/", filenameList)
