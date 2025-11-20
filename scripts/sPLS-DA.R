# ============================================================================
# sPLS-DA ANALYSIS WITH LOADING PLOTS, SCORE PLOTS, AND STATISTICAL TESTS
# ============================================================================
# Author: David Rovelli
# Date: 2025-11-03
# Purpose: Complete sPLS-DA analysis including visualizations and statistics
# ============================================================================

# Install packages (uncomment if needed)
# if (!requireNamespace("BiocManager", quietly = TRUE))    
#   install.packages("BiocManager")
# BiocManager::install("mixOmics")
# install.packages(c("ggplot2", "dplyr", "tidyr", "svglite", "ggpubr", "rstatix"))

# ============================================================================
# LOAD LIBRARIES
# ============================================================================
library(mixOmics)
library(ggplot2)
library(dplyr)
library(tidyr)
library(svglite)
library(ggpubr)
library(rstatix)

# ============================================================================
# LOAD DATA
# ============================================================================
hot_infuses <- read.csv("C:\\Users\\david\\OneDrive - Università degli Studi di Parma\\PhD\\Projects\\Coffee_leaf_infuses\\data\\processed\\others\\hot_scaled_metaboanalyst.csv", row.names = 1)
cold_infuses <- read.csv("C:\\Users\\david\\OneDrive - Università degli Studi di Parma\\PhD\\Projects\\Coffee_leaf_infuses\\data\\processed\\others\\cold_scaled_metaboanalyst.csv", row.names = 1)

# ============================================================================
# FUNCTION: PREPARE DATA AND CREATE sPLS-DA MODEL
# ============================================================================
splsda_model <- function(data, keepX = c(15, 15)) {
  X <- t(data)
  
  # Create the response factor
  sample_names <- rownames(X)
  Y <- sapply(strsplit(sample_names, "_"), `[`, 1)
  Y <- as.factor(Y)
  
  # sPLS-DA model
  splsda_result <- splsda(X, Y, keepX = keepX)
  
  return(splsda_result)
}

# ============================================================================
# FUNCTION: CREATE LOADING PLOTS
# ============================================================================
plot_loadings_points <- function(splsda_obj, comp = 1, ndisplay = 25, 
                                 title_text = "Loadings", outcome_colors = NULL) {
  
  # Extract loadings for the component
  loadings_data <- splsda_obj$loadings$X[, comp]
  
  # Extract X matrix and Y factor
  X_data <- splsda_obj$X
  Y_data <- splsda_obj$Y
  
  # Create dataframe with loadings
  loadings_df <- data.frame(
    variable = names(loadings_data),
    loading = as.numeric(loadings_data)
  )
  
  # Sort by absolute loading and keep top ndisplay
  loadings_df <- loadings_df %>%
    arrange(desc(abs(loading))) %>%
    slice(1:ndisplay)
  
  # For each variable, find which outcome has the maximum mean loading contribution
  outcome_assignment <- sapply(loadings_df$variable, function(var) {
    var_idx <- which(colnames(X_data) == var)
    if(length(var_idx) > 0) {
      # Calculate mean value for each outcome
      mean_values <- tapply(X_data[, var_idx], Y_data, mean, na.rm = TRUE)
      # Return the outcome with maximum absolute value
      names(mean_values)[which.max(abs(mean_values))]
    } else {
      return(NA)
    }
  })
  
  loadings_df$outcome <- outcome_assignment
  
  # Order variables by loading value for better visualization
  loadings_df <- loadings_df %>%
    mutate(variable = factor(variable, levels = rev(variable)))
  
  # Define default colors for outcomes (professional palette)
  if(is.null(outcome_colors)) {
    outcomes_unique <- unique(Y_data)
    outcome_colors <- c(
      "#E75480", "#0B7285", "#2E7D32", "#F57C00", "#1565C0"
    )
    names(outcome_colors) <- levels(outcomes_unique)
  }
  
  # Create the plot
  p <- ggplot(loadings_df, aes(x = loading, y = variable, fill = outcome)) +
    geom_point(size = 4, stroke = 1.5, color = "black", shape = 21) +
    geom_segment(aes(xend = 0, yend = variable, color = NA), size = 0.3, color = "#CCCCCC", linetype = "dotted") +
    scale_fill_manual(
      values = outcome_colors,
      name = "Legend",
      na.value = "#CCCCCC"
    ) +
    theme_minimal() +
    labs(
      title = title_text,
      subtitle = paste("Top", ndisplay, "variables - Component", comp),
      x = paste("Loadings", comp),
      y = NULL
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#555555", margin = margin(b = 15)),
      axis.title.x = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid.major.x = element_line(color = "#E0E0E0", size = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  return(p)
}

# ============================================================================
# FUNCTION: CREATE SCORE PLOTS
# ============================================================================
plot_scoreplot <- function(splsda_obj, comp_x = 1, comp_y = 2,
                          title_text = "Score Plot", outcome_colors = NULL) {
  
  # Extract scores
  scores <- splsda_obj$variates$X
  Y_data <- splsda_obj$Y
  
  # Create dataframe with scores
  scores_df <- data.frame(
    comp_x = scores[, comp_x],
    comp_y = scores[, comp_y],
    outcome = Y_data
  )
  
  # Define default colors for outcomes
  if(is.null(outcome_colors)) {
    outcomes_unique <- unique(Y_data)
    outcome_colors <- c(
      "#E75480", "#0B7285", "#2E7D32", "#F57C00", "#1565C0"
    )
    names(outcome_colors) <- levels(outcomes_unique)
  }
  
  # Calculate explained variance
  explained_var_x <- round(splsda_obj$prop_expl_var$X[comp_x] * 100, 2)
  explained_var_y <- round(splsda_obj$prop_expl_var$X[comp_y] * 100, 2)
  
  # Create the plot
  p <- ggplot(scores_df, aes(x = comp_x, y = comp_y, color = outcome, fill = outcome)) +
    geom_point(size = 4, stroke = 1.5, alpha = 0.8, shape = 21) +
    scale_color_manual(
      values = outcome_colors,
      name = "Legend"
    ) +
    scale_fill_manual(
      values = outcome_colors,
      name = "Legend"
    ) +
    theme_minimal() +
    labs(
      title = title_text,
      subtitle = paste("Component", comp_x, "vs Component", comp_y),
      x = paste("Component", comp_x, "-", explained_var_x, "% Explained Variance"),
      y = paste("Component", comp_y, "-", explained_var_y, "% Explained Variance")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#CCCCCC", size = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#CCCCCC", size = 0.5) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#555555", margin = margin(b = 15)),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "#E0E0E0", size = 0.2),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(15, 15, 15, 15),
      aspect.ratio = 1
    )
  
  return(p)
}

# ============================================================================
# FUNCTION: GENERATE AND SAVE PLS-DA PLOTS
# ============================================================================
generate_and_save_plots <- function(splsda_obj, dataset_name, outcome_colors, 
                                   output_dir = "scripts/outputs/sPLS-DA") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Generating loading plots...\n")
  
  # Generate loading plots for both components
  plot_comp1 <- plot_loadings_points(splsda_obj, comp = 1, ndisplay = 25,
                                     title_text = "Loadings 1",
                                     outcome_colors = outcome_colors)
  
  plot_comp2 <- plot_loadings_points(splsda_obj, comp = 2, ndisplay = 25,
                                     title_text = "Loadings 2",
                                     outcome_colors = outcome_colors)
  
  # Save loading plots in both PNG and SVG formats
  for (comp in 1:2) {
    plot <- if(comp == 1) plot_comp1 else plot_comp2
    
    # PNG
    png_path <- file.path(output_dir, paste0(dataset_name, "_loadings_comp", comp, "_points.png"))
    ggsave(png_path, plot, width = 10, height = 8, dpi = 300)
    
    # SVG
    svg_path <- file.path(output_dir, paste0(dataset_name, "_loadings_comp", comp, "_points.svg"))
    ggsave(svg_path, plot, width = 10, height = 8, dpi = 300)
  }

  cat("Generating score plots...\n")

  # Generate scoreplot
  scoreplot_comp1_2 <- plot_scoreplot(splsda_obj, comp_x = 1, comp_y = 2,
                                      title_text = "Score Plot (Comp 1 vs 2)",
                                      outcome_colors = outcome_colors)
  
  # Save scoreplot in both PNG and SVG formats
  # PNG
  png_path_score <- file.path(output_dir, paste0(dataset_name, "_scoreplot_comp1_comp2.png"))
  ggsave(png_path_score, scoreplot_comp1_2, width = 10, height = 8, dpi = 300)
  
  # SVG
  svg_path_score <- file.path(output_dir, paste0(dataset_name, "_scoreplot_comp1_comp2.svg"))
  ggsave(svg_path_score, scoreplot_comp1_2, width = 10, height = 8, dpi = 300)
  
  cat("Displaying plots...\n")
  
  # Display plots
  print(plot_comp1)
  print(plot_comp2)
  print(scoreplot_comp1_2)
  
  cat("\n✓ Plots saved for:", dataset_name, "\n")
  cat("  -", dataset_name, "_loadings_comp1/2_points.png/svg\n")
  cat("  -", dataset_name, "_scoreplot_comp1_comp2.png/svg\n\n")
}

# ============================================================================
# FUNCTION: CREATE BOXPLOTS FOR SELECTED MOLECULES
# ============================================================================
plot_molecules_boxplot <- function(data, molecule_list, dataset_name, 
                                   outcome_colors = NULL, nrow = 2, ncol = 3) {
  
  # Transpose data for easier access
  X <- t(data)
  
  # Create the response factor
  sample_names <- rownames(X)
  Y <- sapply(strsplit(sample_names, "_"), `[`, 1)
  Y <- as.factor(Y)
  
  # Create dataframe for plotting
  plot_data <- data.frame()
  
  for (molecule in molecule_list) {
    if (molecule %in% colnames(X)) {
      temp_df <- data.frame(
        molecule = molecule,
        value = X[, molecule],
        outcome = Y
      )
      plot_data <- rbind(plot_data, temp_df)
    } else {
      cat("Warning: Molecule", molecule, "not found in data\n")
    }
  }
  
  # Define default colors if not provided
  if(is.null(outcome_colors)) {
    outcomes_unique <- unique(Y)
    outcome_colors <- c(
      "#E75480", "#0B7285", "#2E7D32", "#F57C00", "#1565C0"
    )
    names(outcome_colors) <- levels(outcomes_unique)
  }
  
  # Create faceted boxplot
  p <- ggplot(plot_data, aes(x = outcome, y = value, fill = outcome)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2, 
                  outlier.stroke = 1, color = "black") +
    geom_jitter(width = 0.2, size = 2, alpha = 0.5, shape = 21, 
                color = "black", stroke = 0.5) +
    scale_fill_manual(
      values = outcome_colors,
      name = "Legend"
    ) +
    facet_wrap(~molecule, scales = "free_y", nrow = nrow, ncol = ncol) +
    theme_minimal() +
    labs(
      title = paste("Selected Key Odorant Compounds -", dataset_name),
      x = "Samples",
      y = "Scaled Intensity"
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "#E0E0E0", size = 0.2),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  return(p)
}

# ============================================================================
# FUNCTION: SAVE MOLECULES PLOTS
# ============================================================================
save_molecules_plots <- function(plot_obj, dataset_name, 
                                molecule_names = NULL,
                                suffix = "",
                                output_dir = "scripts/outputs/sPLS-DA",
                                width = 14, height = 10) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate filename
  filename_base <- paste0(dataset_name, "_molecules_boxplot", suffix)
  if (!is.null(molecule_names) && length(molecule_names) > 0) {
    short_names <- substr(molecule_names[1:min(3, length(molecule_names))], 1, 10)
    filename_base <- paste0(filename_base, "_", paste(short_names, collapse = "_"))
  }
  
  # PNG
  png_path <- file.path(output_dir, paste0(filename_base, ".png"))
  ggsave(png_path, plot_obj, width = width, height = height, dpi = 300)
  cat("✓ PNG saved:", basename(png_path), "\n")
  
  # SVG
  svg_path <- file.path(output_dir, paste0(filename_base, ".svg"))
  ggsave(svg_path, plot_obj, width = width, height = height, dpi = 300)
  cat("✓ SVG saved:", basename(svg_path), "\n")
}

# ============================================================================
# DEFINE COLORS (CONSISTENT FOR ALL ANALYSES)
# ============================================================================
hot_outcome_colors <- c(
  "H.M.EC" = "#E75480",
  "H.M.T5" = "#0B7285",
  "H.MA" = "#2E7D32",
  "H.N.GE.WS" = "#F57C00",
  "H.WI" = "#1565C0"
)

cold_outcome_colors <- c(
  "C.M.EC" = "#E75480",
  "C.M.T5" = "#0B7285",
  "C.MA" = "#2E7D32",
  "C.N.GE.WS" = "#F57C00",
  "C.WI" = "#1565C0"
)

# ============================================================================
# DEFINE MOLECULES OF INTEREST
# ============================================================================
hot_molecules_of_interest <- c(
  "Benzeneethanol", 
  "Pentanal", 
  "Linalool", 
  "Linalool oxide", 
  "trans-beta-Ionone", 
  "Benzene methanol", 
  "beta-Cyclocitral",
  "Octanoic acid", 
  "3,5-Octadien-2-one", 
  "Hexanal", 
  "2-Hexenal, (E)-", 
  "Benzaldehyde", 
  "Nonanal", 
  "Methyl salicylate"
)


cold_molecules_of_interest <- c(
  "Benzeneethanol",
  "Linalool",
  "D-Limonene", 
  "trans-beta-Ionone",
  "Benzene methanol",
  "beta-Cyclocitral",
  "3,5-Octadien-2-one",
  "2-Hexenal, (E)-",
  "Benzaldehyde",
  "Nonanal", 
  "Methyl salicylate"
)

# ============================================================================
# RUN HOT INFUSES - sPLS-DA ANALYSIS
# ============================================================================
cat("\n", strrep("=", 80), "\n")
cat("PROCESSING HOT INFUSES - sPLS-DA\n")
cat(strrep("=", 80), "\n\n")

splsda_hot <- splsda_model(hot_infuses, keepX = c(15, 15))
generate_and_save_plots(splsda_hot, "hot", hot_outcome_colors)

# ============================================================================
# RUN COLD INFUSES - sPLS-DA ANALYSIS
# ============================================================================
cat("\n", strrep("=", 80), "\n")
cat("PROCESSING COLD INFUSES - sPLS-DA\n")
cat(strrep("=", 80), "\n\n")

splsda_cold <- splsda_model(cold_infuses, keepX = c(15, 15))
generate_and_save_plots(splsda_cold, "cold", cold_outcome_colors)

# ============================================================================
# BOXPLOTS FOR HOT INFUSES
# ============================================================================
cat("\n", strrep("=", 80), "\n")
cat("GENERATING BOXPLOTS FOR HOT INFUSES\n")
cat(strrep("=", 80), "\n\n")

hot_molecules_plot <- plot_molecules_boxplot(hot_infuses, 
                                             hot_molecules_of_interest,
                                             "HOT",
                                             hot_outcome_colors,
                                             nrow = 3, ncol = 5)
print(hot_molecules_plot)
save_molecules_plots(hot_molecules_plot, "hot", hot_molecules_of_interest, width = 18, height = 12)

# ============================================================================
# BOXPLOTS FOR COLD INFUSES
# ============================================================================
cat("\n", strrep("=", 80), "\n")
cat("GENERATING BOXPLOTS FOR COLD INFUSES\n")
cat(strrep("=", 80), "\n\n")

cold_molecules_plot <- plot_molecules_boxplot(cold_infuses, 
                                              cold_molecules_of_interest,
                                              "COLD",
                                              cold_outcome_colors,
                                              nrow = 3, ncol = 5)
print(cold_molecules_plot)
save_molecules_plots(cold_molecules_plot, "cold", cold_molecules_of_interest, width = 18, height = 12)


# ============================================================================
# ANALYSIS SUMMARY
# ============================================================================
cat("\n", strrep("=", 80), "\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat(strrep("=", 80), "\n\n")

cat("OUTPUT GENERATED:\n")
cat("  ✓ Molecules boxplots without p-values (hot & cold)\n")

cat("Files saved in: scripts/outputs/sPLS-DA/\n")
cat(strrep("=", 80), "\n")