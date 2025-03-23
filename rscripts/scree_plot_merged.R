library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("M:/T1DGC/USERS/dam8mt/data/MEGA")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eigenval = read.delim("PCA/mega_pca_b37_pruned_unrelated_KGref.eigenval", header = F)

# Convert column into rows
eigenval = as.data.frame(t(eigenval))

# Name the columns
colnames(eigenval) = c(paste0("PC", 1:10))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data frame
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Sum of eigenvalues (PC1:PC10)
sum = rowSums(eigenval)

# Proportion of variance - divide each eigenvalue (PC1:PC10) by sum of eigenvalues
var_explained = mapply('/', eigenval, sum)

# A data frame containing the PCs and the variance explained by each PC
var_explained_df = data.frame(PC=(1:10), var_explained=var_explained)

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

dataList = list(var_explained_df)
titleList = list("Scree plot merged")
filenameList = list("scree_plot_merged.png")

mapply(screePlot, dataList, titleList, "Graphs/PCA_merged/", filenameList)
