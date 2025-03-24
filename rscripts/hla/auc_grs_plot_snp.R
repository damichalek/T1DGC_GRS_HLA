library(ggplot2)
library(pROC)
library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("M:/T1DGC/USERS/dam8mt/data/MEGA")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define functions to read and label files
# %s is a placeholder for strings.
# %d is a placeholder for integers.
# %f is a placeholder for floating-point numbers.
read_and_label = function(anc1, anc2) {
  path = sprintf("HLA/GRS/output_AMR/chr6_%s_SNP_grs_in_%sgrs.txt", anc1, anc2)
  data = read.table(path, header = T)
  assign(sprintf("%sin%s", anc1, anc2), data, envir = .GlobalEnv)
}

# Load data
for (anc1 in c("AFR", "AMR", "EUR", "FIN")) {
  for (anc2 in c("AFR", "AMR", "EUR", "FIN")) {
    read_and_label(anc1, anc2)
  }
}

read_and_label_merged = function(anc2) {
  path = sprintf("HLA/GRS_merged/output/chr6_merged_SNP_grs_in_%sgrs.txt", anc2)
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

# Function to add GRS_label and Population to the data
add_labels = function(anc1, anc2) {
  data_name = sprintf("%sin%s", anc1, anc2)
  data = get(data_name)
  
  # Add the GRS_label and Population columns
  data = data %>% 
    mutate(GRS_label = anc1,
           Population = anc2)
  
  # Assign the modified data back to the original variable
  assign(data_name, data, envir = .GlobalEnv)
}

# Loop through all combinations of ancestries to add the labels
for (anc1 in c("AFR", "AMR", "EUR", "FIN")) {
  for (anc2 in c("AFR", "AMR", "EUR", "FIN")) {
    add_labels(anc1, anc2)
  }
}

# Function to add GRS_label and Population to the merged data
add_labels_merged = function(anc2) {
  data_name = sprintf("ALLin%s", anc2)
  data = get(data_name)
  
  # Add the GRS_label and Population columns
  data = data %>% 
    mutate(GRS_label = "ALL",
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
  AFRinAFR, AFRinAMR, AFRinEUR, AFRinFIN,
  AMRinAFR, AMRinAMR, AMRinEUR, AMRinFIN,
  EURinAFR, EURinAMR, EURinEUR, EURinFIN,
  FINinAFR, FINinAMR, FINinEUR, FINinFIN,
  ALLinAFR, ALLinAMR, ALLinEUR, ALLinFIN
)

# Combine all datasets into one using bind_rows
combined = bind_rows(data_list)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AUC plot (each ancestry separately)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aucPlot = function(data, pred, obs, grs_label, population, path) {
  
  # Define colors for each population
  population_colors <- c(
    AFR = "#85D4E3",
    AMR = "#F4B5BD",
    EUR = "#CDC08C",
    FIN = "#FAD77B"
  )
  
  # Filter data for the specific GRS and Population
  filtered_data = data %>% 
    filter(GRS_label == grs_label) %>% 
    filter(Population == population)
  
  # Calculate ROC and AUC using pROC package
  roc_obj = roc(filtered_data[[obs]], filtered_data[[pred]])
  auc_val = auc(roc_obj)
  
  # Create a data frame for plotting ROC curve
  roc_data = data.frame(
    tpr = roc_obj$sensitivities,   # True positive rate (Sensitivity)
    fpr = 1 - roc_obj$specificities # False positive rate (1 - Specificity)
  )
  
  # Plot using ggplot2
  plot = ggplot(roc_data, aes(x = fpr, y = tpr)) +
    geom_line(color = population_colors[population], linewidth = 1.5) +  # Use specific color
    geom_abline(slope = 1, intercept = 0, linetype = "solid") +          # Diagonal line
    ggtitle(paste("T1D Prediction:", grs_label, "GRS in", population, "ancestry")) +  # Title
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_bw() +                                     # Minimal theme
    xlab("False Positive Rate (1 - Specificity)") +                         # X-axis label
    ylab("True Positive Rate (Sensitivity)") +                          # Y-axis label
    annotate("text", x = 0.8, y = 0.1,                    # Add AUC value
             label = paste("AUC =", round(auc_val, 3))) +
    annotate("text", x = 0.8, y = 0.05,                   # Number of positive cases
             label = paste("Npres =", sum(filtered_data[[obs]] == 1))) +
    annotate("text", x = 0.8, y = 0,                      # Number of negative cases
             label = paste("Nabs =", sum(filtered_data[[obs]] == 0))) +
    theme(
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  
  # Save the plot
  filename = sprintf("%sauc_%s_grs_in_%s_population_snp.png", path, grs_label, population)
  ggsave(filename, plot, width = 15, height = 15, units = "cm", dpi = 300, scaling = 1.1)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate and save plots (each ancestry separately)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Unique GRS labels and populations
grs_labels = unique(combined$GRS_label)
populations = unique(combined$Population)

# Loop through each combination of GRS and Population
for (grs_label in grs_labels) {
  for (population in populations) {
    # Skip invalid combinations (e.g., "ALL" GRS with "merged" population doesn't exist)
    if (grs_label == "ALL" && !(population %in% c("AFR", "AMR", "EUR", "FIN"))) next
    
    # Generate the plot
    aucPlot(
      data = combined, 
      pred = "ScaledGRS", 
      obs = "Status", 
      grs_label = grs_label, 
      population = population, 
      path = "HLA/Graphs/AUC/SNP/"
    )
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AUC plot (combined by ancestry)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define a function to plot AUC for multiple GRS labels in each ancestry population
aucPlot_by_ancestry = function(data, pred, obs, path) {
  
  # Define the GRS labels you want to plot
  grs_labels = c("AFR", "AMR", "EUR", "FIN", "ALL")
  
  # Define colors for each GRS label
  grs_colors <- c(
    AFR = "#85D4E3",
    AMR = "#CDC08C",
    EUR = "#F4B5BD",
    FIN = "#C2C2F0",
    ALL = "#FAD77B"
  )
  
  # Define populations for which we want to generate the plots
  populations = c("AFR", "AMR", "EUR", "FIN")
  
  # Loop through each population (AFR, AMR, EUR, FIN)
  for (population in populations) {
    
    # Create an empty list to store ROC data for each GRS label
    roc_data_list = list()
    
    # Loop through each GRS label and calculate ROC and AUC
    for (grs_label in grs_labels) {
      
      # Filter data for the specific GRS and population
      filtered_data = data %>% 
        filter(GRS_label == grs_label) %>% 
        filter(Population == population)
      
      # Calculate ROC using pROC package
      roc_obj = roc(filtered_data[[obs]], filtered_data[[pred]])
      
      # Create a data frame for plotting ROC curve with the GRS label and color
      roc_data = data.frame(
        tpr = roc_obj$sensitivities,   # True positive rate (Sensitivity)
        fpr = 1 - roc_obj$specificities, # False positive rate (1 - Specificity)
        GRS_label = grs_label           # Add the GRS label for color and legend
      )
      
      # Append the ROC data for this GRS label to the list
      roc_data_list[[grs_label]] = roc_data
    }
    
    # Combine all ROC data frames into one
    combined_roc_data = bind_rows(roc_data_list)
    
    # Define the GRS labels in the desired order
    grs_labels = c("AFR", "AMR", "EUR", "FIN", "ALL")
    
    # Ensure factor levels are correctly ordered
    combined_roc_data$GRS_label = factor(combined_roc_data$GRS_label, levels = grs_labels)
    
    # Plot using ggplot2
    plot = ggplot(combined_roc_data, aes(x = fpr, y = tpr, color = GRS_label)) +
      geom_line(linewidth = 1.5) +  # Plot all ROC curves with different colors
      geom_abline(slope = 1, intercept = 0, linetype = "solid") +  # Diagonal line
      #ggtitle(paste("T1D Prediction: AUC for Different GRS Labels in", population, "Ancestry")) +  # Title
      #theme(plot.title = element_text(hjust = 0.5)) +
      theme_bw() +  # Minimal theme
      xlab("False Positive Rate (1 - Specificity)") +  # X-axis label
      ylab("True Positive Rate (Sensitivity)") +  # Y-axis label
      scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2)) +  # Set x-axis limits and remove padding
      scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2)) +  # Set y-axis limits and remove padding
      scale_color_manual(values = grs_colors) +  # Custom colors for each GRS label
      labs(color = "GRS model") +  # Add legend title
      theme(
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),  # Make legend title bold
        legend.position = "bottom"  # Place legend at the bottom
      )
    
    # Save the plot for the current population
    filename = sprintf("%sauc_combined_grs_in_%s_population_snp.png", path, population)
    ggsave(filename, plot, width = 15, height = 15, units = "cm", dpi = 300, scaling = 1.15)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate and save plots for all ancestries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aucPlot_by_ancestry(
  data = combined, 
  pred = "ScaledGRS", 
  obs = "Status", 
  path = "HLA/Graphs/AUC/SNP/"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AUC plot (combined by ancestry with facet_wrap)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define a function to plot AUC for multiple GRS labels in all ancestries
aucPlot_combined_ancestries = function(data, pred, obs, path) {
  
  # Define the GRS labels you want to plot for each ancestry
  grs_labels = c("AFR", "AMR", "EUR", "FIN", "ALL")
  
  # Create an empty list to store ROC data for each GRS label and ancestry
  roc_data_list = list()
  
  # Loop through each ancestry and GRS label to calculate ROC and AUC
  for (population in c("AFR", "AMR", "EUR", "FIN")) {
    
    # Loop through each GRS label
    for (grs_label in grs_labels) {
      
      # Filter data for the specific GRS label and population
      filtered_data = data %>% 
        filter(GRS_label == grs_label) %>% 
        filter(Population == population)
      
      # Skip if there's not enough data to calculate ROC
      if (nrow(filtered_data) < 2) {
        warning(paste("Not enough data for ROC calculation for", grs_label, "in", population, "population"))
        next
      }
      
      # Calculate ROC using pROC package
      roc_obj = tryCatch({
        roc(filtered_data[[obs]], filtered_data[[pred]])
      }, error = function(e) {
        warning(paste("Error in ROC calculation for", grs_label, "in", population, "population:", e$message))
        return(NULL)
      })
      
      # Skip if ROC calculation failed
      if (is.null(roc_obj)) next
      
      # Create a data frame for plotting ROC curve with the GRS label, population, and color
      roc_data = data.frame(
        tpr = roc_obj$sensitivities,   # True positive rate (Sensitivity)
        fpr = 1 - roc_obj$specificities, # False positive rate (1 - Specificity)
        GRS_label = grs_label,           # Add the GRS label for color and legend
        Population = population         # Add the population for faceting
      )
      
      # Append the ROC data for this GRS label and population to the list
      roc_data_list[[paste(grs_label, population)]] = roc_data
    }
  }
  
  # Combine all ROC data frames into one
  combined_roc_data = bind_rows(roc_data_list)
  
  # Define the GRS labels in the desired order
  grs_labels = c("AFR", "AMR", "EUR", "FIN", "ALL")
  
  # Ensure factor levels are correctly ordered
  combined_roc_data$GRS_label = factor(combined_roc_data$GRS_label, levels = grs_labels)
  
  # Define colors for each GRS label
  grs_colors <- c(
    AFR = "#85D4E3",
    AMR = "#CDC08C",
    EUR = "#F4B5BD",
    FIN = "#C2C2F0",
    ALL = "#FAD77B"
  )
  
  # Plot using ggplot2
  plot = ggplot(combined_roc_data, aes(x = fpr, y = tpr, color = GRS_label)) +
    geom_line(linewidth = 0.8) +  # Plot all ROC curves with different colors
    geom_abline(slope = 1, intercept = 0, linetype = "solid") +  # Diagonal line
    theme_bw() +  # Minimal theme
    xlab("False Positive Rate (1 - Specificity)") +  # X-axis label
    ylab("True Positive Rate (Sensitivity)") +  # Y-axis label
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2)) +  # Set x-axis limits and remove padding
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2)) +  # Set y-axis limits and remove padding
    scale_color_manual(values = grs_colors) +  # Custom colors for each GRS label
    labs(color = bquote(bold("T1D GRS"["HLA-SNP"] ~ "models"))) +  # Add legend title
    theme(
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),  # Make legend title bold
      legend.position = "bottom",  # Place legend at the bottom
      strip.background = element_rect(fill = "#273046", color = "#273046"),  # Change facet_wrap background color
      strip.text = element_text(color = "white", face = "bold")) +  # Change facet text to white for contrast
    facet_wrap(~ Population, ncol = 2)  # Create separate facets for each population
  
  # Save the combined plot
  filename = sprintf("%sauc_combined_grs_in_all_populations_snp.png", path)
  ggsave(filename, plot, width = 20, height = 20, units = "cm", dpi = 300, scaling = 1.15)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate and save plot (combined by ancestry)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aucPlot_combined_ancestries(
  data = combined, 
  pred = "ScaledGRS", 
  obs = "Status", 
  path = "HLA/Graphs/AUC/SNP/"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AUC plot (combined by ancestry with facet_wrap and AUC values)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define a function to plot AUC for multiple GRS labels in all ancestries
aucPlot_combined_ancestries_auc = function(data, pred, obs, path) {
  
  # Define the GRS labels you want to plot for each ancestry
  grs_labels = c("AFR", "AMR", "EUR", "FIN", "ALL")
  
  # Create an empty list to store ROC data for each GRS label and ancestry
  roc_data_list = list()
  
  # Create a list to store AUC values
  auc_values = list()
  
  # Loop through each ancestry and GRS label to calculate ROC and AUC
  for (population in c("AFR", "AMR", "EUR", "FIN")) {
    
    # Loop through each GRS label
    for (grs_label in grs_labels) {
      
      # Filter data for the specific GRS label and population
      filtered_data = data %>% 
        filter(GRS_label == grs_label) %>% 
        filter(Population == population)
      
      # Skip if there's not enough data to calculate ROC
      if (nrow(filtered_data) < 2) {
        warning(paste("Not enough data for ROC calculation for", grs_label, "in", population, "population"))
        next
      }
      
      # Calculate ROC using pROC package
      roc_obj = tryCatch({
        roc(filtered_data[[obs]], filtered_data[[pred]])
      }, error = function(e) {
        warning(paste("Error in ROC calculation for", grs_label, "in", population, "population:", e$message))
        return(NULL)
      })
      
      # Skip if ROC calculation failed
      if (is.null(roc_obj)) next
      
      # Store AUC value
      auc_value = roc_obj$auc[[1]]
      auc_values[[paste(grs_label, population)]] = auc_value
      
      # Create a data frame for plotting ROC curve with the GRS label, population, and color
      roc_data = data.frame(
        tpr = roc_obj$sensitivities,   # True positive rate (Sensitivity)
        fpr = 1 - roc_obj$specificities, # False positive rate (1 - Specificity)
        GRS_label = grs_label,           # Add the GRS label for color and legend
        Population = population         # Add the population for faceting
      )
      
      # Append the ROC data for this GRS label and population to the list
      roc_data_list[[paste(grs_label, population)]] = roc_data
    }
  }
  
  # Combine all ROC data frames into one
  combined_roc_data = bind_rows(roc_data_list)
  
  # Define the GRS labels in the desired order
  grs_labels = c("AFR", "AMR", "EUR", "FIN", "ALL")
  
  # Ensure factor levels are correctly ordered
  combined_roc_data$GRS_label = factor(combined_roc_data$GRS_label, levels = grs_labels)
  
  # Define colors for each GRS label
  grs_colors <- c(
    AFR = "#85D4E3",
    AMR = "#CDC08C",
    EUR = "#F4B5BD",
    FIN = "#C2C2F0",
    ALL = "#FAD77B"
  )
  
  # Plot using ggplot2
  plot = ggplot(combined_roc_data, aes(x = fpr, y = tpr, color = GRS_label)) +
    geom_line(linewidth = 0.8) +  # Plot all ROC curves with different colors
    geom_abline(slope = 1, intercept = 0, linetype = "solid") +  # Diagonal line
    theme_bw() +  # Minimal theme
    xlab("False Positive Rate (1 - Specificity)") +  # X-axis label
    ylab("True Positive Rate (Sensitivity)") +  # Y-axis label
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2)) +  # Set x-axis limits and remove padding
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2)) +  # Set y-axis limits and remove padding
    scale_color_manual(values = grs_colors) +  # Custom colors for each GRS label
    labs(color = "GRS model") +  # Add legend title
    theme(
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),  # Make legend title bold
      legend.position = "bottom",  # Place legend at the bottom
      strip.background = element_rect(fill = "#273046", color = "#273046"),  # Change facet_wrap background color
      strip.text = element_text(color = "white", face = "bold")) +  # Change facet text to white for contrast
    facet_wrap(~ Population, ncol = 2) +  # Create separate facets for each population
    geom_text(data = combined_roc_data %>%
                mutate(GRS_label = factor(GRS_label, levels = c("FIN", "EUR", "AMR", "AFR", "ALL"))) %>%  # Correct factor ordering
                group_by(Population, GRS_label) %>%
                summarise(auc_label = paste0("AUC~GRS[", GRS_label, "] ==", round(auc_values[[paste(GRS_label[1], Population[1])]], 2))),
              aes(x = 0.6, y = 0.08 + 0.05 * as.numeric(factor(GRS_label)), label = auc_label), size = 3, color = "black", hjust = 0, # Adjust position for each GRS
              parse = TRUE)  # Enable parsing of mathematical expressions
  
  # Save the combined plot
  filename = sprintf("%sauc_combined_grs_in_all_populations_auc_snp.png", path)
  ggsave(filename, plot, width = 20, height = 20, units = "cm", dpi = 300, scaling = 1.15)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate and save plot (combined by ancestry)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aucPlot_combined_ancestries_auc(
  data = combined, 
  pred = "ScaledGRS", 
  obs = "Status", 
  path = "HLA/Graphs/AUC/SNP/"
)
