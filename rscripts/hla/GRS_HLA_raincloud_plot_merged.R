library(tidyverse)
library(ghibli)
library(ggdist)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("M:/T1DGC/USERS/dam8mt/data/MEGA")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

read_and_label_merged = function(anc2) {
  path = sprintf("HLA/GRS_merged/output/chr6_merged_HLA_grs_in_%sgrs.txt", anc2)
  data = read.table(path, header = T)
  assign(sprintf("ALLin%s", anc2), data, envir = .GlobalEnv)
}

# Load merged data
for (anc2 in c("AFR", "AMR", "EUR", "FIN")) {
  read_and_label_merged(anc2)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add GRS and Population variable
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to add GRS_label and Population to the merged data
add_labels_merged = function(anc2) {
  data_name = sprintf("ALLin%s", anc2)
  data = get(data_name)
  
  # Add the GRS_label and Population columns
  data = data %>% 
    mutate(GRS_label = "GRS ALL",
           Population = anc2)
  
  # Assign the modified data back to the original variable
  assign(data_name, data, envir = .GlobalEnv)
}

# Add labels to merged data
for (anc2 in c("AFR", "AMR", "EUR", "FIN")) {
  add_labels_merged(anc2)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create a list of all datasets
data_list = list(
  ALLinAFR, ALLinAMR, ALLinEUR, ALLinFIN
)

# Combine all datasets into one using bind_rows
combined = bind_rows(data_list)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRS raincloud plots in different populations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to plot GRS raincloud plot
raincloudPlot = function(data, path, filename, ancestry, ext){
  data %>% 
    filter(GRS_label == paste0("GRS ", ancestry)) %>%  # Filter for specific GRS label
    mutate(Status = factor(ifelse(Status == 0, "Control", "Case"), levels = c("Control", "Case"))) %>%  # Set Control first
    ggplot(aes(x = 1, y = GRS, fill = Status)) +  # Use constant x for overlap
    stat_halfeye(justification = -0.2, .width = 0, point_colour = NA, alpha = 0.6) +
    geom_boxplot(width = 0.12, outlier.color = NA, alpha = 0.6) +
    stat_dots(side = "left", justification = 1.1, binwidth = NA, aes(color = Status), alpha = 0.6) +
    #geom_half_point(side = "l", range_scale = 0, shape = 95, size = 1.5, alpha = .3) +
    scale_fill_manual(values = c("#F5CDB4", "#9A8822")) +
    scale_color_manual(values = c("#F5CDB4", "#9A8822")) +
    theme_bw() +
    xlab("") +
    ylab(bquote(bold(T1D~GRS)[bold(HLA-Allele-.(ancestry))])) + # Dynamic Subscript with bquote()
    #ylab(paste0("T1D ", ancestry, "-derived GRS")) +  # Dynamic Y-axis label without "GRS" +
    #ylim(-10, 5) +
    theme(
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      axis.text.y = element_blank(),  # Remove x-axis text
      axis.ticks.y = element_blank(),  # Remove x-axis ticks
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#273046", color = "#273046"),  # Change facet_wrap background color
      strip.text = element_text(color = "white", face = "bold")) +  # Change facet text to white for contrast
    facet_wrap(~ Population) +  # Facet by Population to get separate plots for each ancestry
    coord_flip()
  
  # Save the plot
  # Dynamically update the filename based on GRS label
  ggsave(paste0(path, filename, "_", ancestry, ext), width = 18, height = 13, units = "cm", dpi = 300, scaling = 0.9)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ancestries = c("ALL")  # List of ancestries
dataList = list(combined)
mapply(raincloudPlot, dataList, "HLA/Graphs/GRS/", "grs_raincloud_plot_hla", ancestries, ".png")
