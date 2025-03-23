library(ggplot2)
library(ggpubr) #ggarrange
library(dplyr)
library(RColorBrewer)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Choose color for plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# View a single RColorBrewer palette by specifying its name
display.brewer.pal(n = 12, name = 'Paired')

# Hexadecimal color specification 
brewer.pal(n = 12, name = "Paired")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("M:/T1DGC/USERS/dam8mt/data/MEGA")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load PCA results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pca = read.delim("PCA/PCs/mega_b37_pcs_AMR.txt", header = T)
pca$AFF = as.factor(pca$AFF)
pca = pca %>% 
  filter(AFF !=-9)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Scatter plot of PCA results (colored by affection status)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# COLOR SELECTION
# Select 3rd and 4th color of a twelve-color Paired palette:
myColors <- brewer.pal(12, "Paired")[c(8,10)]

plotPCA = function(data, x, y, labx, laby, path, filename){
  data %>% 
    ggplot(aes(x = x, y = y, color = AFF)) + 
    geom_point(size = 2.5) +
    labs(x = labx, y = laby) +
    theme_bw() +
    scale_color_manual(values = myColors)
  
  ggsave(paste0(path,filename), width = 17, height = 15, units="cm", dpi=300)
}

# Save plots

plotPCA(pca, pca$PC1, pca$PC1, "PC1", "PC1", "Graphs/PCA/aff/AMR/", "pc1_pc1.png")
plotPCA(pca, pca$PC1, pca$PC2, "PC1", "PC2", "Graphs/PCA/aff/AMR/", "pc1_pc2.png")
plotPCA(pca, pca$PC1, pca$PC3, "PC1", "PC3", "Graphs/PCA/aff/AMR/", "pc1_pc3.png")
plotPCA(pca, pca$PC1, pca$PC4, "PC1", "PC4", "Graphs/PCA/aff/AMR/", "pc1_pc4.png")
plotPCA(pca, pca$PC1, pca$PC5, "PC1", "PC5", "Graphs/PCA/aff/AMR/", "pc1_pc5.png")

plotPCA(pca, pca$PC2, pca$PC1, "PC2", "PC1", "Graphs/PCA/aff/AMR/", "pc2_pc1.png")
plotPCA(pca, pca$PC2, pca$PC2, "PC2", "PC2", "Graphs/PCA/aff/AMR/", "pc2_pc2.png")
plotPCA(pca, pca$PC2, pca$PC3, "PC2", "PC3", "Graphs/PCA/aff/AMR/", "pc2_pc3.png")
plotPCA(pca, pca$PC2, pca$PC4, "PC2", "PC4", "Graphs/PCA/aff/AMR/", "pc2_pc4.png")
plotPCA(pca, pca$PC2, pca$PC5, "PC2", "PC5", "Graphs/PCA/aff/AMR/", "pc2_pc5.png")

plotPCA(pca, pca$PC3, pca$PC1, "PC3", "PC1", "Graphs/PCA/aff/AMR/", "pc3_pc1.png")
plotPCA(pca, pca$PC3, pca$PC2, "PC3", "PC2", "Graphs/PCA/aff/AMR/", "pc3_pc2.png")
plotPCA(pca, pca$PC3, pca$PC3, "PC3", "PC3", "Graphs/PCA/aff/AMR/", "pc3_pc3.png")
plotPCA(pca, pca$PC3, pca$PC4, "PC3", "PC4", "Graphs/PCA/aff/AMR/", "pc3_pc4.png")
plotPCA(pca, pca$PC3, pca$PC5, "PC3", "PC5", "Graphs/PCA/aff/AMR/", "pc3_pc5.png")

plotPCA(pca, pca$PC4, pca$PC1, "PC4", "PC1", "Graphs/PCA/aff/AMR/", "pc4_pc1.png")
plotPCA(pca, pca$PC4, pca$PC2, "PC4", "PC2", "Graphs/PCA/aff/AMR/", "pc4_pc2.png")
plotPCA(pca, pca$PC4, pca$PC3, "PC4", "PC3", "Graphs/PCA/aff/AMR/", "pc4_pc3.png")
plotPCA(pca, pca$PC4, pca$PC4, "PC4", "PC4", "Graphs/PCA/aff/AMR/", "pc4_pc4.png")
plotPCA(pca, pca$PC4, pca$PC5, "PC4", "PC5", "Graphs/PCA/aff/AMR/", "pc4_pc5.png")

plotPCA(pca, pca$PC5, pca$PC1, "PC5", "PC1", "Graphs/PCA/aff/AMR/", "pc5_pc1.png")
plotPCA(pca, pca$PC5, pca$PC2, "PC5", "PC2", "Graphs/PCA/aff/AMR/", "pc5_pc2.png")
plotPCA(pca, pca$PC5, pca$PC3, "PC5", "PC3", "Graphs/PCA/aff/AMR/", "pc5_pc3.png")
plotPCA(pca, pca$PC5, pca$PC4, "PC5", "PC4", "Graphs/PCA/aff/AMR/", "pc5_pc4.png")
plotPCA(pca, pca$PC5, pca$PC5, "PC5", "PC5", "Graphs/PCA/aff/AMR/", "pc5_pc5.png")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine plots into one
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SavePlotAsVariable = function(data, x, y, labx, laby){
  data %>% 
    ggplot(aes(x = x, y = y, color = AFF)) + 
    geom_point(size = 2.5) +
    labs(x = labx, y = laby) +
    theme_bw() +
    scale_color_manual(values = myColors)
}

# Save plots as variables

p11 = SavePlotAsVariable(pca, pca$PC1, pca$PC1, "PC1", "PC1")
p12 = SavePlotAsVariable(pca, pca$PC1, pca$PC2, "PC1", "PC2")
p13 = SavePlotAsVariable(pca, pca$PC1, pca$PC3, "PC1", "PC3")
p14 = SavePlotAsVariable(pca, pca$PC1, pca$PC4, "PC1", "PC4")
p15 = SavePlotAsVariable(pca, pca$PC1, pca$PC5, "PC1", "PC5")

p21 = SavePlotAsVariable(pca, pca$PC2, pca$PC1, "PC2", "PC1")
p22 = SavePlotAsVariable(pca, pca$PC2, pca$PC2, "PC2", "PC2")
p23 = SavePlotAsVariable(pca, pca$PC2, pca$PC3, "PC2", "PC3")
p24 = SavePlotAsVariable(pca, pca$PC2, pca$PC4, "PC2", "PC4")
p25 = SavePlotAsVariable(pca, pca$PC2, pca$PC5, "PC2", "PC5")

p31 = SavePlotAsVariable(pca, pca$PC3, pca$PC1, "PC3", "PC1")
p32 = SavePlotAsVariable(pca, pca$PC3, pca$PC2, "PC3", "PC2")
p33 = SavePlotAsVariable(pca, pca$PC3, pca$PC3, "PC3", "PC3")
p34 = SavePlotAsVariable(pca, pca$PC3, pca$PC4, "PC3", "PC4")
p35 = SavePlotAsVariable(pca, pca$PC3, pca$PC5, "PC3", "PC5")

p41 = SavePlotAsVariable(pca, pca$PC4, pca$PC1, "PC4", "PC1")
p42 = SavePlotAsVariable(pca, pca$PC4, pca$PC2, "PC4", "PC2")
p43 = SavePlotAsVariable(pca, pca$PC4, pca$PC3, "PC4", "PC3")
p44 = SavePlotAsVariable(pca, pca$PC4, pca$PC4, "PC4", "PC4")
p45 = SavePlotAsVariable(pca, pca$PC4, pca$PC5, "PC4", "PC5")

p51 = SavePlotAsVariable(pca, pca$PC5, pca$PC1, "PC5", "PC1")
p52 = SavePlotAsVariable(pca, pca$PC5, pca$PC2, "PC5", "PC2")
p53 = SavePlotAsVariable(pca, pca$PC5, pca$PC3, "PC5", "PC3")
p54 = SavePlotAsVariable(pca, pca$PC5, pca$PC4, "PC5", "PC4")
p55 = SavePlotAsVariable(pca, pca$PC5, pca$PC5, "PC5", "PC5")

p = ggarrange(p11, p12, p13, p14, p15,
              p21, p22, p23, p24, p25,
              p31, p32, p33, p34, p35,
              p41, p42, p43, p44, p45,
              p51, p52, p53, p54, p55,
              ncol = 5, nrow = 5)

ggsave("Graphs/PCA/aff/AMR/all_pc_aff_AMR.png", plot = p, width = 85, height = 75, units="cm", dpi=300)
