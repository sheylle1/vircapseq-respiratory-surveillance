library(tidyverse)
library(readxl)
library(pROC)
library(ggplot2)


INPUT_DIR <- "data/processed"
OUTPUT_DIR <- "figures"
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
  
  pcr_data <- read_excel(file.path(input_dir, pcr_file)) %>%
    select(sample_id, Pathogen, PCR_Result) %>%                 # Keep only needed columns
    filter(!is.na(PCR_Result), PCR_Result != -1)                # Remove NA and "Not tested"
  
  continuous_data <- read_csv(file.path(input_dir, continuous_file), show_col_types = FALSE) %>%
    select(sample_id, Pathogen, coverage_percentage, rpm, depth_of_coverage, Contam) %>%
    mutate( across(c(coverage_percentage, rpm, depth_of_coverage), as.numeric), Contam = ifelse(is.na(Contam), "NO", toupper(Contam)))
  
  roc_data <- inner_join(pcr_data, continuous_data, by = c("sample_id", "Pathogen")) %>%
    drop_na(PCR_Result, coverage_percentage, rpm, depth_of_coverage)
  
  roc_data <- roc_data %>%
    mutate(
      coverage_percentage = ifelse(Contam == "YES", 0, coverage_percentage),
      rpm = ifelse(Contam == "YES", 0, rpm),
      depth_of_coverage = ifelse(Contam == "YES", 0, depth_of_coverage)
    )
  
  # Optional: Show dataset breakdown if original_dataset exists
  if ("original_dataset" %in% colnames(pcr_data)) {
    dataset_summary <- roc_data %>%
      left_join(pcr_data %>% select(sample_id, Pathogen, original_dataset) %>% distinct(),
                by = c("sample_id", "Pathogen")) %>%
      group_by(original_dataset) %>%
      summarise(
        total = n(),
        positive = sum(PCR_Result == 1),
        negative = sum(PCR_Result == 0)
      )
    print(dataset_summary)  # Add print to show summary
  } 
  
  return(roc_data)  # Make sure to return the data
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
  
  auc_values <- tibble( Metric = c("Genome coverage", "RPM"), AUC = c(auc_coverage, auc_rpm))
  
  # Order
  roc_list <- list(RPM = roc_rpm, `Genome coverage` = roc_coverage)
  
  list(
    roc_objects = roc_list,
    auc_values = auc_values,
    auc_coverage = auc_coverage,
    auc_rpm = auc_rpm,
    auc_depth = auc_depth
  )
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
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 13),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 18, hjust = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.25)
    )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  base_name <- file.path(output_dir, "Figure1_A")
  
  ggsave(paste0(base_name, ".png"), plot = p, 
         width = 10, height = 7, dpi = 300, bg = "white")
  ggsave(paste0(base_name, ".pdf"), plot = p, 
         width = 10, height = 7, device = cairo_pdf, bg = "white")
  
  return(invisible(list(plot = p, auc = roc_results$auc_values)))
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
  
  roc_data <- load_roc_data(INPUT_DIR, PCR_FILE, CONTINUOUS_FILE)
  
  # Check if we have both classes
  if (sum(roc_data$PCR_Result == 1) == 0) {
    stop("No positive samples found in the data!")
  }
  if (sum(roc_data$PCR_Result == 0) == 0) {
    stop("No negative samples found in the data!")
  }
  
  roc_results <- generate_roc_curves(roc_data) 
  create_roc_plot(roc_results, OUTPUT_DIR) 
}

# Run
main()