# # Install BiocManager if not installed 
# if (!requireNamespace("BiocManager", quietly = TRUE))    
# install.packages("BiocManager")
# # Install mixOmics
# BiocManager::install("mixOmics")
# # Install other required CRAN packages
# install.packages(c("ggplot2", "dplyr", "tidyr", "svglite"))

# Load libraries
library(mixOmics)
library(ggplot2)
library(dplyr)
library(tidyr)
library(svglite)

# Load data
hot_infuses <- read.csv("C:\\Users\\david\\OneDrive - Università degli Studi di Parma\\PhD\\Projects\\Coffee_leaf_infuses\\data\\processed\\others\\hot_scaled_metaboanalyst.csv", row.names = 1)
cold_infuses <- read.csv("C:\\Users\\david\\OneDrive - Università degli Studi di Parma\\PhD\\Projects\\Coffee_leaf_infuses\\data\\processed\\others\\cold_scaled_metaboanalyst.csv", row.names = 1)

# Function to prepare data and create sPLS-DA model
splsda_model <- function(data, keepX = c(50, 30)) {
  X <- t(data)
  
  # Create the response factor
  sample_names <- rownames(X)
  Y <- sapply(strsplit(sample_names, "_"), `[`, 1)
  Y <- as.factor(Y)
  
  # sPLS-DA model
  splsda_result <- splsda(X, Y, keepX = keepX)
  
  return(splsda_result)
}

# Function to create loading plots
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
      name = "Outcome",
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

# Function to create scoreplot
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
      name = "Outcome"
    ) +
    scale_fill_manual(
      values = outcome_colors,
      name = "Outcome"
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

# Function to generate and save plots (loadings + scoreplots)
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
  cat("  - hot_loadings_comp1/2_points.png/svg\n")
  cat("  - hot_scoreplot_comp1_comp2.png/svg\n\n")
}

# Define colors (consistent for all analyses)
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
# RUN HOT INFUSES
# ============================================================================
cat("Processing HOT infuses...\n")
splsda_hot <- splsda_model(hot_infuses, keepX = c(50, 30))
generate_and_save_plots(splsda_hot, "hot", hot_outcome_colors)

# ============================================================================
# RUN COLD INFUSES
# ============================================================================
cat("Processing COLD infuses...\n")
splsda_cold <- splsda_model(cold_infuses, keepX = c(50, 30))
generate_and_save_plots(splsda_cold, "cold", cold_outcome_colors)

cat("Analyses completed!\n")