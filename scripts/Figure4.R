#!/usr/bin/env Rscript

library(tidyverse)
library(readxl)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(scales)

# --- Helper functions --------------------------------------------------

convert_detected <- function(x) {
  dplyr::case_when(
    stringr::str_to_lower(x) == "detected" ~ 1,
    stringr::str_to_lower(x) == "not detected" ~ 0,
    stringr::str_to_lower(x) == "not tested" ~ -1,
    TRUE ~ as.numeric(x)
  )
}

prepare_main_data <- function(main_file) {
  read_csv(main_file, show_col_types = FALSE) %>%
    mutate(
      Contam = ifelse(is.na(Contam) | Contam == "", "NO", toupper(Contam)),
      across(c(coverage_percentage, rpm, depth_of_coverage, ct_values), as.numeric)
    )
}

load_metadata <- function(metadata_file) {
  metadata <- read_csv(metadata_file)
  
  if ("fileName" %in% colnames(metadata)) {
    metadata <- metadata %>% rename(sample_id = fileName)
  }
  
  metadata <- metadata %>%
    mutate(
      totalReads_clean = as.numeric(gsub(",", "", totalReads)),
      total_reads_category = case_when(
        totalReads_clean < 100000 ~ "<100K",
        totalReads_clean >= 100000 & totalReads_clean < 500000 ~ "100K-500K",
        totalReads_clean >= 500000 & totalReads_clean < 1000000 ~ "500K-1M",
        totalReads_clean >= 1000000 & totalReads_clean < 2000000 ~ "1M-2M",
        totalReads_clean >= 2000000 & totalReads_clean < 5000000 ~ "2M-5M",
        totalReads_clean >= 5000000 & totalReads_clean < 10000000 ~ "5M-10M",
        totalReads_clean >= 10000000 ~ ">10M",
        TRUE ~ NA_character_
      )
    )
  
  read_categories <- c("<100K", "100K-500K", "500K-1M", "1M-2M", "2M-5M", "5M-10M", ">10M")
  
  list(metadata = metadata, read_categories = read_categories)
}

create_ct_figure <- function(main_data, metadata_info, output_dir) {
  
  all_data <- main_data %>%
    left_join(metadata_info$metadata, by = "sample_id") %>%
    filter(!is.na(ct_values) & ct_values > 0) %>%
    mutate(
      Ct_bin = cut(ct_values, 
                   breaks = c(-Inf, 20, 25, 30, 35, Inf),
                   labels = c("Ct <20", "Ct 20-25", "Ct 25-30", "Ct 30-35", "Ct >35")),
      total_reads_category = factor(total_reads_category, 
                                    levels = metadata_info$read_categories)
    ) %>%
    filter(!is.na(total_reads_category)) %>%
    mutate(
      depth_of_coverage = ifelse(depth_of_coverage == 0 & Contam != "YES", 0.1, depth_of_coverage),
      coverage_percentage = ifelse(coverage_percentage == 0 & Contam != "YES", 0.1, coverage_percentage),
      is_zero_point = (depth_of_coverage == 0.1 & coverage_percentage == 0.1),
      is_contam_point = (Contam == "YES")
    )
  
  max_depth_global <- max(all_data$depth_of_coverage[all_data$depth_of_coverage > 0.1], na.rm = TRUE)
  xmax_box_global <- max_depth_global * 1.35
  
  if (max_depth_global > 100000) {
    axis_breaks <- c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000)
    axis_labels <- c("0", "1", "10", "100", "1K", "10K", "100K", "1M")
  } else if (max_depth_global > 10000) {
    axis_breaks <- c(0.1, 1, 10, 100, 1000, 10000, 100000)
    axis_labels <- c("0", "1", "10", "100", "1K", "10K", "100K")
  } else if (max_depth_global > 1000) {
    axis_breaks <- c(0.1, 1, 10, 100, 1000, 10000)
    axis_labels <- c("0", "1", "10", "100", "1K", "10K")
  } else {
    axis_breaks <- c(0.1, 1, 10, 100, 1000)
    axis_labels <- c("0", "1", "10", "100", "1K")
  }
  
  red_box_data <- tibble(
    xmin = 10,
    xmax = xmax_box_global,
    ymin = 80,
    ymax = 105
  )
  
  all_data <- all_data %>%
    mutate(Pathogen = ifelse(Pathogen == "Enterovirus / Rhinovirus (HRV)",
                             "Rhinovirus (HRV)", Pathogen))
  
  pathogens <- sort(unique(all_data$Pathogen))
  
  bold_palette <- c(
    "#E2C044", "#cd7cac", "#72b043", "#F46036",
    "#D7263D", "#662E9B", "#0D3B66", "#2EC4B6"
  )
  
  pathogen_colors <- bold_palette[1:length(pathogens)]
  names(pathogen_colors) <- pathogens
  
  zero_point_counts <- all_data %>%
    filter(is_zero_point) %>%
    group_by(total_reads_category, Ct_bin) %>%
    summarise(zero_point_count = n(), .groups = 'drop')
  
  all_data <- all_data %>%
    left_join(zero_point_counts, by = c("total_reads_category", "Ct_bin")) %>%
    mutate(
      zero_point_count = replace_na(zero_point_count, 0),
      point_color = ifelse(is_zero_point & Contam != "YES" & zero_point_count >= 2, 
                           "gray50", Pathogen)
    )
  
  row_counts <- all_data %>%
    group_by(total_reads_category) %>%
    summarise(row_n = n(), .groups = 'drop')
  
  col_counts <- all_data %>%
    group_by(Ct_bin) %>%
    summarise(col_n = n(), .groups = 'drop')
  
  all_data <- all_data %>%
    left_join(row_counts, by = "total_reads_category") %>%
    left_join(col_counts, by = "Ct_bin") %>%
    mutate(
      total_reads_category_label = paste0(total_reads_category, "\n(n=", row_n, ")"),
      Ct_bin_label = paste0(Ct_bin, "\n(n=", col_n, ")")
    )
  
  all_data$total_reads_category_label <- factor(
    all_data$total_reads_category_label,
    levels = unique(all_data$total_reads_category_label[order(all_data$total_reads_category)])
  )
  
  all_data$Ct_bin_label <- factor(
    all_data$Ct_bin_label,
    levels = unique(all_data$Ct_bin_label[order(all_data$Ct_bin)])
  )
  
  point_counts <- all_data %>%
    filter(is_zero_point) %>%
    group_by(total_reads_category_label, Ct_bin_label, depth_of_coverage, coverage_percentage) %>%
    summarise(count = n(), .groups = 'drop') %>%
    filter(count >= 1)
  
  p <- ggplot(all_data, aes(x = depth_of_coverage, y = coverage_percentage)) +
    geom_rect(data = red_box_data,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
              color = "red", fill = NA, linewidth = 0.3,
              inherit.aes = FALSE) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_point(data = subset(all_data, !is_contam_point),
               aes(color = point_color), alpha = 0.7, size = 1.2) +
    geom_point(data = subset(all_data, is_contam_point),
               aes(fill = point_color),
               color = "red", size = 2.0, stroke = 0.45,
               shape = 21, alpha = 0.9) +
    geom_text(data = point_counts,
              aes(x = 0.1, y = 0.1, label = paste0("n=", count)),
              color = "black", size = 3.2, 
              vjust = -0.35, hjust = 0.4,
              inherit.aes = FALSE,
              check_overlap = FALSE) +
    facet_grid(total_reads_category_label ~ Ct_bin_label, 
               scales = "fixed", space = "fixed") +
    scale_fill_manual(values = c(pathogen_colors, "gray50" = "gray50"), 
                      name = "Pathogen Detected",
                      breaks = names(pathogen_colors),
                      na.value = "grey50",
                      guide = "none") +
    scale_color_manual(values = c(pathogen_colors, "gray50" = "gray50"), 
                       name = "Pathogen Detected",
                       breaks = names(pathogen_colors),
                       na.value = "grey50") +
    scale_y_continuous(labels = percent_format(scale = 1), limits = c(0, 105)) +
    scale_x_log10(breaks = axis_breaks, labels = axis_labels, limits = c(0.09, xmax_box_global)) +
    labs(x = "Depth of coverage", y = "Genome coverage (%)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = 12, face = "bold", vjust = -0.4),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 8, face = "bold"),
      strip.text.y = element_text(size = 8, face = "bold"),
      strip.background = element_rect(color = "black", fill = "grey95", linewidth = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.spacing = unit(0.02, "lines"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      legend.box.margin = margin(t = 15, b = 10),
      legend.background = element_rect(fill = "white", color = "white", linewidth = 0.5),
      plot.margin = margin(20, 10, 20, 10)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))
  
  ggsave(filename = file.path(output_dir, "Figure4.png"), plot = p,
         width = 9, height = 11, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(filename = file.path(output_dir, "Figure4.pdf"), plot = p,
         width = 9, height = 11, device = cairo_pdf, bg = "white", limitsize = FALSE)
  
  return(list(plot = p, data = all_data))
}

output_dir <- "figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

metadata_file <- "data/processed/sequencing_reads.csv"
if (file.exists(metadata_file)) {
  metadata_info <- load_metadata(metadata_file)
} else {
  stop("Metadata file not found: ", metadata_file)
}

main_file <- "data/processed/meta_data.csv"
if (file.exists(main_file)) {
  main_data <- prepare_main_data(main_file)
} else {
  stop("Main data file not found: ", main_file)
}

if (nrow(main_data) > 0) {
  results <- create_ct_figure(main_data, metadata_info, output_dir)
}