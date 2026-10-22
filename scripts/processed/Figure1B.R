#!/usr/bin/env Rscript
library(tidyverse)
library(readxl)
library(patchwork)
library(cowplot)
library(scales)
library(grid)
library(gtable)
library(svglite)

INPUT_DIR <- "data/processed/"
OUTPUT_DIR <- "results/Figures"
PCR_FILE <- "data_PCR.xlsx"
CONTINUOUS_FILE <- "data_detections.csv"

COVERAGE_THRESHOLDS <- c(1, 10, 20, 50, 80, 90)
RPM_THRESHOLDS <- c(1, 10, 100, 1000)

COLOR_LOW  <- "#D73027"
COLOR_MID  <- "#FEE08B"
COLOR_HIGHLOW1 <- "#D9F0D3"
COLOR_HIGHLOW2 <- "#A6DBA0"
COLOR_HIGHMID1 <- "#8CBC5B"
COLOR_HIGHMID2 <- "#6DAF38"
COLOR_HIGHMID3 <- "#4F9C78"
COLOR_HIGHMID4 <- "#0A8456"
COLOR_HIGH  <- "#006D3A"

HEADER_FILL <- "#f2f2f2"
ROWLBL_FILL <- "#f2f2f2"
HEADER_BORDER <- "black"
INNER_BORDER  <- "white"

METRIC_SCALES <- list(
  Sensitivity = list(min = 0.57, mid = 0.80, max = 1.00),
  Specificity = list(min = 0.57, mid = 0.8,  max = 1.00),
  Accuracy  = list(min = 0.57, mid = 0.8,  max = 1.00)
)

METRICS <- c("Sensitivity", "Specificity", "Accuracy")

COLOR_VALS <- c(0.57, 0.8, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 1.0)
ALL_COLORS <- c(COLOR_LOW, COLOR_MID, COLOR_HIGHLOW1, COLOR_HIGHLOW2,
                COLOR_HIGHMID1, COLOR_HIGHMID2, COLOR_HIGHMID3,
                COLOR_HIGHMID4, COLOR_HIGH)

# ── Data loading ──────────────────────────────────────────────────────────────
load_data <- function(input_dir, pcr_file, continuous_file) {
  # Read PCR data
  pcr_data <- read_excel(file.path(input_dir, pcr_file)) %>%
    select(sample_id, Pathogen, PCR_Result) %>%
    filter(!is.na(PCR_Result), PCR_Result != -1)  # Remove NA and -1
  
  # Read continuous data
  continuous_data <- read_csv(file.path(input_dir, continuous_file), show_col_types = FALSE) %>%
    select(sample_id, Pathogen, coverage_percentage, rpm, Contam) %>%
    mutate(
      across(c(coverage_percentage, rpm), as.numeric),
      Contam = ifelse(is.na(Contam), "NO", toupper(Contam))
    )
  
  # Get all sample-pathogen combinations
  all_combos <- pcr_data %>% select(sample_id, Pathogen) %>% distinct()
  
  # Complete the continuous data with all combinations
  continuous_complete <- all_combos %>%
    left_join(continuous_data, by = c("sample_id", "Pathogen")) %>%
    mutate(
      rpm = replace_na(rpm, 0),
      coverage_percentage = replace_na(coverage_percentage, 0),
      Contam = replace_na(Contam, "NO")
    )
  
  # Join with PCR data
  merged_data <- inner_join(pcr_data, continuous_complete,  by = c("sample_id", "Pathogen"))
  
  # Print summary info
  cat("Total samples:", nrow(merged_data), "\n")
  cat("PCR Positive:", sum(merged_data$PCR_Result == 1), "\n")
  cat("PCR Negative:", sum(merged_data$PCR_Result == 0), "\n")
  cat("Contaminated samples:", sum(merged_data$Contam == "YES"), "\n")
  
  return(merged_data)
}

calculate_metrics <- function(actual, predicted) {
  # Remove NA values
  v  <- !is.na(actual) & !is.na(predicted)
  a  <- actual[v]
  p  <- predicted[v]
  
  TP <- sum(a == 1 & p == 1)
  TN <- sum(a == 0 & p == 0)
  FP <- sum(a == 0 & p == 1)
  FN <- sum(a == 1 & p == 0)
  
  list(
    TP = TP, 
    TN = TN, 
    FP = FP, 
    FN = FN,
    Sensitivity = ifelse((TP+FN) > 0, round(TP/(TP+FN), 3), NA_real_),
    Specificity = ifelse((TN+FP) > 0, round(TN/(TN+FP), 3), NA_real_),
    Accuracy    = ifelse((TP+TN+FP+FN) > 0,  round((TP+TN)/(TP+TN+FP+FN), 3), NA_real_)
  )
}

compute_grid <- function(merged_data, thresholds, column, type_label) {
  map_dfr(thresholds, function(thresh) {
    # For contaminated samples, set prediction to 0 (negative)
    # For non-contaminated, use threshold
    pred <- ifelse(merged_data$Contam == "YES", 0, ifelse(merged_data[[column]] >= thresh, 1, 0))
    
    m <- calculate_metrics(merged_data$PCR_Result, pred)
    
    tibble(Threshold   = as.character(thresh),
           TP = m$TP, 
           TN = m$TN, 
           FP = m$FP, 
           FN = m$FN,
           Sensitivity = m$Sensitivity,
           Specificity = m$Specificity,
           Accuracy    = m$Accuracy,
           Metric_Type = type_label)
  })
}

label_fn <- function(x) {
  # Only show values that are not exactly 1.000 or 0.000
  # And not NA
  ifelse(is.na(x) | x == 1.000 | x == 0.000, "", sprintf("%.3f", x))
}

make_col_header <- function(txt) {
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = HEADER_FILL, colour = HEADER_BORDER, linewidth = 0.5) +
    annotate("text", x = 0.5, y = 0.5,  label = txt, size = 3.2, fontface = "bold", family = "sans", colour = "grey10") +
    theme_void() +
    theme(plot.margin = margin(6, 0, 6, 0))  # t/b margin controls header height
}

make_row_label <- function(label) {
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = ROWLBL_FILL, colour = HEADER_BORDER, linewidth = 0.5) +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 3.2, fontface = "bold",= family = "sans", angle = 90,= hjust = 0.5, vjust = 0.5, colour = "grey15") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_tile <- function(data, metric, x_label, show_x_axis,
                      is_coverage = FALSE, is_rpm = FALSE) {
  sc  <- METRIC_SCALES[[metric]]
  lim  <- c(sc$min, sc$max)
  norm_vals <- scales::rescale(COLOR_VALS, from = lim)
  
  # Margin: add gap on right of coverage and left of RPM to create divider space
  tile_margin <- if (is_coverage) margin(0, 4, 0, 0)
                 else if (is_rpm) margin(0, 0, 0, 4)
                 else margin(0, 0, 0, 0)
  
  ggplot(data,
         aes(x = factor(Threshold, levels = unique(Threshold)),  y = 0.5, fill = .data[[metric]])) +
    geom_tile(colour = INNER_BORDER, linewidth = 0.4) +
    geom_text(aes(label = label_fn(.data[[metric]])), colour = "grey15", size = 2.8, fontface = "bold", family = "sans") +
    scale_fill_gradientn(
      colors = ALL_COLORS,
      values  = norm_vals,
      limits = lim,
      na.value = "grey92",
      guide  = "none"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = if (show_x_axis) x_label else NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x  = if (show_x_axis) 
                        element_text(size = 8, margin = margin(t = 2))
                     else element_blank(),
      axis.ticks.x = if (show_x_axis) 
                        element_line(linewidth = 0.3)
                     else element_blank(),
      axis.title.x = element_text(size = 8, face = "bold",  margin = margin(t = 4)),
      panel.grid   = element_blank(),
      panel.border = element_rect(colour = HEADER_BORDER,fill = NA, linewidth = 0.5),
      plot.margin  = tile_margin
    )
}

make_unified_legend <- function() {
  lim  <- c(0.57, 1.0)
  norm_vals <- scales::rescale(COLOR_VALS, from = lim)
  y_steps <- seq(lim[1], lim[2], length.out = 200)
  tile_h  <- (lim[2] - lim[1]) / 199
  
  dummy <- ggplot(tibble(x = 1, y = y_steps, v = y_steps), aes(x = x, y = y, fill = v)) +
    geom_tile(width = 1, height = tile_h, linewidth = 0) +
    scale_fill_gradientn(
      colors = ALL_COLORS,
      values = norm_vals,
      limits = lim,
      breaks = c(0.6, 0.7, 0.8, 0.9, 1.0),
      labels = function(b) sprintf("%.1f", b),
      name = NULL,
      guide = guide_colorbar(
        barwidth = unit(0.5, "cm"),
        barheight  = unit(5,   "cm"),
        title.position = "top",
        label.theme = element_text(size = 8),
        ticks.colour = "grey40",
        frame.colour  = HEADER_BORDER,
        frame.linewidth = 0.5
      )
    ) +
    theme_void() + theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0))
  
  get_legend(dummy)
}

create_figure <- function(coverage_grid, rpm_grid, output_dir) {
  
  # Top header row
  header_gap <- ggplot() + theme_void()
  header_cov <- make_col_header("Genome Coverage (%)")
  header_rpm <- make_col_header("Reads per million (RPM)")
  
  header_row <- header_gap + header_cov + header_rpm +
    plot_layout(nrow = 1, widths = c(0.8, 4, 3.2))
  
  # Data rows
  rows <- vector("list", length(METRICS))
  for (i in seq_along(METRICS)) {
    m         <- METRICS[i]
    is_bottom <- (i == length(METRICS))
    
    tile_cov <- make_tile(coverage_grid, m,  x_label     = "", show_x_axis = is_bottom, is_coverage = TRUE)
    tile_rpm <- make_tile(rpm_grid, m, x_label     = "", show_x_axis = is_bottom, is_rpm      = TRUE)
    row_lbl  <- make_row_label(m)
    
    rows[[i]] <- row_lbl + tile_cov + tile_rpm +
      plot_layout(nrow = 1, widths = c(0.8, 4, 3.2))
  }
  
  body <- wrap_plots(
    header_row,
    rows[[1]], rows[[2]], rows[[3]],
    ncol    = 1,
    heights = c(0.15, 0.28, 0.28, 0.28)
  ) &
  theme(plot.margin = margin(0, 0, 0, 0))
  
  legend_grob <- make_unified_legend()
  
  final <- (body | wrap_elements(legend_grob)) +
    plot_layout(widths = c(9, 0.7)) +
    plot_annotation(
      theme = theme(
        plot.margin = margin(4, 4, 4, 4),
        plot.background = element_rect(fill    = "white", colour  = HEADER_BORDER,blinewidth = 0.5)
      )
    )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  base <- file.path(output_dir, "Figure1_B")
  
  ggsave(paste0(base, ".png"), plot = final,
         width = 170, height = 120, units = "mm", dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), plot = final,
         width = 170, height = 120, units = "mm",
         device = cairo_pdf, bg = "white")
  ggsave(paste0(base, ".svg"), plot = final,
         width = 170, height = 120, units = "mm",
         device = svglite::svglite, bg = "white")
  
  cat("\nFigure saved to:", output_dir, "\n")
  invisible(final)
}

main <- function() {
  # Check if input directory exists, if not use current directory
  if (!dir.exists(INPUT_DIR)) {
    cat("Warning: Directory", INPUT_DIR, "not found. Using current directory.\n")
    INPUT_DIR <- "."
  }
  
  merged_data <- load_data(INPUT_DIR, PCR_FILE, CONTINUOUS_FILE)
  coverage_grid <- compute_grid(merged_data, COVERAGE_THRESHOLDS, "coverage_percentage", "Coverage")
  rpm_grid <- compute_grid(merged_data, RPM_THRESHOLDS, "rpm", "RPM")
  
  # Print the calculated metrics for verification
  cat("\n Coverage Metrics \n")
  print(coverage_grid %>% select(Threshold, Sensitivity, Specificity, Accuracy))
  
  cat("\n RPM Metrics \n")
  print(rpm_grid %>% select(Threshold, Sensitivity, Specificity, Accuracy))
  
  create_figure(coverage_grid, rpm_grid, OUTPUT_DIR)
}

main()