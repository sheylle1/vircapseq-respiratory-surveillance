library(tidyverse)
library(readxl)
library(pROC)
library(ggplot2)
library(svglite)

INPUT_DIR <- "data/processed/"
OUTPUT_DIR <- "results/Figures"
PCR_FILE <- "data_PCR.xlsx"                # Combined PCR file
CONTINUOUS_FILE <- "data_detections.csv"   # Combined detections file

COLOR_BLUE <- "#4575B4"   
COLOR_RED  <- "#D73027" 


convert_detected <- function(x) {
  case_when(
    str_to_lower(x) == "detected" ~ 1,
    str_to_lower(x) == "not detected" ~ 0,
    str_to_lower(x) == "not tested" ~ -1,
    TRUE ~ as.numeric(x)
  )
}


load_roc_data <- function(input_dir, pcr_file, continuous_file) {
  
  # Read PCR data - keep original_dataset for summary
  pcr_data <- read_excel(file.path(input_dir, pcr_file)) %>%
    select(sample_id, Pathogen, PCR_Result, original_dataset) %>%  # Keep original_dataset
    filter(!is.na(PCR_Result))  # Remove NA values
  # Note: PCR_Result should be 0 or 1, not -1 in your data
  
  # Read continuous data
  continuous_data <- read_csv(file.path(input_dir, continuous_file), show_col_types = FALSE) %>%
    select(sample_id, Pathogen, coverage_percentage, rpm, depth_of_coverage, Contam) %>%
    mutate(
      across(c(coverage_percentage, rpm, depth_of_coverage), as.numeric),
      Contam = ifelse(is.na(Contam), "NO", toupper(Contam))
    )
  
  # Join data
  roc_data <- inner_join(pcr_data, continuous_data, by = c("sample_id", "Pathogen")) %>%
    drop_na(PCR_Result, coverage_percentage, rpm, depth_of_coverage)
  
  # Set contaminated samples to 0
  roc_data <- roc_data %>%
    mutate(
      coverage_percentage = ifelse(Contam == "YES", 0, coverage_percentage),
      rpm = ifelse(Contam == "YES", 0, rpm),
      depth_of_coverage = ifelse(Contam == "YES", 0, depth_of_coverage)
    )
  
  # Print dataset summary using the original_dataset column
  if ("original_dataset" %in% colnames(roc_data)) {
    dataset_summary <- roc_data %>%
      group_by(original_dataset) %>%
      summarise(
        total = n(),
        positive = sum(PCR_Result == 1),
        negative = sum(PCR_Result == 0)
      )
    print("Dataset Summary:")
    print(dataset_summary)
  } else {
    # Fallback summary if original_dataset not available
    dataset_summary <- roc_data %>%
      summarise(
        total = n(),
        positive = sum(PCR_Result == 1),
        negative = sum(PCR_Result == 0)
      )
    print("Overall Summary:")
    print(dataset_summary)
  }
  
  return(roc_data)
}

#' Generate ROC curves for all three metrics
generate_roc_curves <- function(roc_data) {
  
  # Validate data has both classes
  if (nrow(roc_data) == 0) { stop("No data available for ROC analysis")}
  if (length(unique(roc_data$PCR_Result)) != 2) {stop("ROC analysis requires both positive (1) and negative (0) PCR results")}
  
  # Compute ROC objects
  roc_coverage <- try(roc(PCR_Result ~ coverage_percentage, data = roc_data), silent = TRUE)
  roc_rpm <- try(roc(PCR_Result ~ rpm, data = roc_data), silent = TRUE)
  roc_depth <- try(roc(PCR_Result ~ depth_of_coverage, data = roc_data), silent = TRUE)
  
  # Extract AUC values
  auc_coverage <- if (!inherits(roc_coverage, "try-error")) auc(roc_coverage) else NA_real_
  auc_rpm <- if (!inherits(roc_rpm, "try-error")) auc(roc_rpm) else NA_real_
  auc_depth <- if (!inherits(roc_depth, "try-error")) auc(roc_depth) else NA_real_
  
  auc_values <- tibble(Metric = c("Genome coverage", "RPM"), AUC = c(auc_coverage, auc_rpm))
  
  # Order
  roc_list <- list(RPM = roc_rpm, `Genome coverage` = roc_coverage)
  
  list(
    roc_objects = roc_list,
    auc_values = auc_values,
    auc_coverage = auc_coverage,
    auc_rpm = auc_rpm,
    auc_depth = auc_depth,
    roc_data = roc_data  # Store original data for threshold calculations
  )
}

#' Calculate confusion matrices at various thresholds
calculate_threshold_tables <- function(roc_data) {
  
  # Define thresholds for coverage
  coverage_thresholds <- c(1, 10, 20, 50, 80, 90)
  
  # Define thresholds for RPM (using reasonable ranges)
  rpm_thresholds <- c(1, 10, 100, 1000)
  
  # Function to calculate confusion matrix at a given threshold
  calculate_confusion <- function(data, metric, threshold) {
    predicted <- ifelse(data[[metric]] >= threshold, 1, 0)
    actual <- data$PCR_Result
    
    TP <- sum(predicted == 1 & actual == 1)
    FP <- sum(predicted == 1 & actual == 0)
    FN <- sum(predicted == 0 & actual == 1)
    TN <- sum(predicted == 0 & actual == 0)
    
    return(tibble(Threshold = as.character(threshold), TP = TP, FP = FP, FN = FN, TN = TN))
  }
  
  # Calculate for coverage thresholds
  coverage_table <- map_dfr(coverage_thresholds, function(thresh) {
    calculate_confusion(roc_data, "coverage_percentage", thresh)
  })
  
  # Calculate for RPM thresholds
  rpm_table <- map_dfr(rpm_thresholds, function(thresh) {
    calculate_confusion(roc_data, "rpm", thresh)
  })
  
  return(list(
    coverage_table = coverage_table,
    rpm_table = rpm_table
  ))
}

create_roc_plot <- function(roc_results, output_dir) {
  
  roc_list <- roc_results$roc_objects
  
  roc_colors <- c(
    "Genome coverage" = COLOR_BLUE, 
    "RPM" = COLOR_RED              
  )
  
  # Legend order from original
  legend_order <- c("Genome coverage", "RPM")
  
  p <- ggroc(roc_list, legacy.axes = TRUE, size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(title = "",
      subtitle = paste(
        "Genome coverage AUC =", round(roc_results$auc_coverage, 3),
        "| RPM AUC =", round(roc_results$auc_rpm, 3)
      ),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    scale_color_manual(
      name = "Metric",
      values = roc_colors,
      breaks = legend_order,
      labels = legend_order
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      plot.subtitle = element_text(size = 10, face = "bold",  hjust = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.25)
    )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  base_name <- file.path(output_dir, "Figure1A")
  
  ggsave(paste0(base_name, ".png"), plot = p, 
         width = 10, height = 7, dpi = 300, bg = "white")
  ggsave(paste0(base_name, ".pdf"), plot = p, 
         width = 10, height = 7, device = cairo_pdf, bg = "white")
  ggsave(paste0(base_name, ".svg"), plot = p, 
       width = 170, height = 120, units = "mm", 
       device = svglite::svglite, bg = "white") 
  
  return(invisible(list(plot = p, auc = roc_results$auc_values)))
}

#' Print threshold tables to console and optionally save to files
print_threshold_tables <- function(threshold_tables, output_dir = NULL) {
  
  cat("\n COVERAGE THRESHOLDS\n")
  print(threshold_tables$coverage_table)
  
  cat("\n RPM THRESHOLDS\n")
  print(threshold_tables$rpm_table)
  
  # Optionally save to CSV files
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    write_csv(threshold_tables$coverage_table, 
              file.path(output_dir, "coverage_thresholds.csv"))
    write_csv(threshold_tables$rpm_table, 
              file.path(output_dir, "rpm_thresholds.csv"))
    
    cat("\nThreshold tables saved to:", output_dir, "\n")
  }
}

#' Calculate and display summary statistics for the data
calculate_summary_stats <- function(roc_data) {
  
  summary_stats <- roc_data %>%
    summarise(
      total_samples = n(),
      positive_pcr = sum(PCR_Result == 1),
      negative_pcr = sum(PCR_Result == 0),
      coverage_mean = mean(coverage_percentage, na.rm = TRUE),
      coverage_median = median(coverage_percentage, na.rm = TRUE),
      coverage_sd = sd(coverage_percentage, na.rm = TRUE),
      rpm_mean = mean(rpm, na.rm = TRUE),
      rpm_median = median(rpm, na.rm = TRUE),
      rpm_sd = sd(rpm, na.rm = TRUE)
    )
  
  print(summary_stats)
  
  return(summary_stats)
}

main <- function() {
  
  # Validate input files exist
  pcr_path <- file.path(INPUT_DIR, PCR_FILE)
  continuous_path <- file.path(INPUT_DIR, CONTINUOUS_FILE)
  
  if (!file.exists(pcr_path)) {
    stop(sprintf("PCR file not found: %s\nExpected path: %s", pcr_path, pcr_path))
  }
  if (!file.exists(continuous_path)) {
    stop(sprintf("Continuous metrics file not found: %s\nExpected path: %s", 
                 continuous_path, continuous_path))
  }
  
  # Load data
  roc_data <- load_roc_data(INPUT_DIR, PCR_FILE, CONTINUOUS_FILE)
  
  # Check if we have both classes
  if (sum(roc_data$PCR_Result == 1) == 0) {
    stop("No positive samples found in the data!")
  }
  if (sum(roc_data$PCR_Result == 0) == 0) {
    stop("No negative samples found in the data!")
  }
  
  cat("\nNumber of samples after filtering:", nrow(roc_data))
  cat("\nPositive samples:", sum(roc_data$PCR_Result == 1))
  cat("\nNegative samples:", sum(roc_data$PCR_Result == 0))
  cat("\n")
  
  # Generate ROC results
  roc_results <- generate_roc_curves(roc_data) 
  
  # Create ROC plot
  create_roc_plot(roc_results, OUTPUT_DIR) 
  
  # Calculate threshold tables
  threshold_tables <- calculate_threshold_tables(roc_data)
  
  # Print threshold tables
  print_threshold_tables(threshold_tables, OUTPUT_DIR)
  
  # Calculate and display summary statistics
  summary_stats <- calculate_summary_stats(roc_data)
  
  # Return everything as a list for potential further use
  return(list(
    roc_results = roc_results,
    threshold_tables = threshold_tables,
    summary_stats = summary_stats,
    roc_data = roc_data
  ))
}

# Run
results <- main()