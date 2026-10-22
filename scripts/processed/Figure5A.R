library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readxl)
library(svglite)

sheets   <- excel_sheets("data/processed/data_coinfect.xlsx")
datasets <- list()
for (i in seq_along(sheets)) {
  datasets[[i]] <- read_excel("data/processed/data_coinfect.xlsx", sheet = sheets[i])
  names(datasets)[i] <- sheets[i]
}

TARGET_SHEETS <- c("Sheet2", "Sheet3", "Sheet4")


AGE_COLOURS <- c(
  "<1 years"    = "#A5C73A",  # green
  "1-5 years"   = "#FFCC00",  # yellow
  "5-18 years"  = "#FF2E2E",  # red
  "18-60 years" = "#33AACC",  # blue
  ">60 years"   = "#BB5FD9",  # purple
  "Unknown"     = "#9E2A5A"   # maroon
)

categorize_age <- function(age) {
  if (is.na(age)) return("Unknown")
  age <- as.numeric(age)
  if      (age < 1)               return("<1 years")
  else if (age >= 1  & age < 5)   return("1-5 years")
  else if (age >= 5  & age < 18)  return("5-18 years")
  else if (age >= 18 & age <= 60) return("18-60 years")
  else if (age > 60)              return(">60 years")
  else                            return("Unknown")
}

create_double_plot <- function(df_long, title, dataset_name = NULL) {

  binary_matrix <- df_long %>%
    mutate(Present = 1) %>%
    pivot_wider(id_cols     = c(sample, Age_Category),
                names_from  = Pathogen,
                values_from = Present,
                values_fill = list(Present = 0)) %>%
    as.data.frame()

  pathogen_cols <- setdiff(colnames(binary_matrix), c("sample", "Age_Category"))
  if (length(pathogen_cols) == 0) return(NULL)

  binary_matrix$CombinationID <- apply(
    binary_matrix[, pathogen_cols, drop = FALSE], 1,
    function(x) paste(x, collapse = "")
  )

  combo_totals <- binary_matrix %>%
    group_by(CombinationID) %>%
    summarise(total_n = n(), .groups = "drop")

  combination_counts <- binary_matrix %>%
    group_by(CombinationID, Age_Category) %>%
    summarise(n = n(), .groups = "drop") %>%
    left_join(
      binary_matrix %>%
        select(-sample, -Age_Category) %>%
        distinct(CombinationID, .keep_all = TRUE),
      by = "CombinationID"
    ) %>%
    rowwise() %>%
    mutate(num_pathogens = sum(c_across(all_of(pathogen_cols)))) %>%
    ungroup() %>%
    filter(num_pathogens >= 2)

  if (nrow(combination_counts) == 0) return(NULL)

  age_levels <- names(AGE_COLOURS)
  age_levels <- age_levels[age_levels %in% unique(combination_counts$Age_Category)]

  combination_counts <- combination_counts %>%
    mutate(Age_Category = factor(Age_Category, levels = age_levels)) %>%
    left_join(combo_totals, by = "CombinationID") %>%
    arrange(Age_Category, desc(total_n), CombinationID)

  # Column order: lowest age group of the combination first, then by total count
  combo_labels <- combination_counts %>%
    distinct(CombinationID, .keep_all = TRUE) %>%
    mutate(id = factor(as.character(seq_len(n())),levels = as.character(seq_len(n())))) %>%
    select(id, CombinationID, total_n, Age_Category)

  combination_counts <- combination_counts %>%
    left_join(combo_labels %>% select(CombinationID, id), by = "CombinationID")

  combo_order <- levels(combo_labels$id)
  n_combos <- length(combo_order)
  combination_counts$id <- factor(combination_counts$id, levels = combo_order)

  pathogens_present <- combination_counts %>%
    select(all_of(pathogen_cols)) %>%
    summarise(across(everything(), ~any(. == 1))) %>%
    pivot_longer(everything(), names_to = "Pathogen", values_to = "present") %>%
    filter(present) %>%
    pull(Pathogen)

  df_long_plot <- combination_counts %>%
    pivot_longer(cols = all_of(pathogen_cols), names_to = "Pathogen", values_to = "Present") %>%
    filter(Pathogen %in% pathogens_present)

  pathogen_order <- sort(unique(pathogens_present))
  if ("CMV" %in% pathogen_order) {
    pathogen_order <- c(pathogen_order[pathogen_order != "CMV"], "CMV")
  }
  df_long_plot$Pathogen <- factor(df_long_plot$Pathogen, levels = pathogen_order)
  df_long_plot$id  <- factor(df_long_plot$id, levels = combo_order)

  # ── Co-infection rate per age group ──────────────────────────────────────
  total_samples_per_age <- binary_matrix %>%
    group_by(Age_Category) %>%
    summarise(total_samples = n(), .groups = "drop")

  coinfected_samples <- binary_matrix %>%
    mutate(num_pathogens = rowSums(select(., all_of(pathogen_cols)))) %>%
    filter(num_pathogens >= 2) %>%
    group_by(Age_Category) %>%
    summarise(coinfected = n(), .groups = "drop")

  age_stats <- total_samples_per_age %>%
    left_join(coinfected_samples, by = "Age_Category") %>%
    mutate(coinfected = ifelse(is.na(coinfected), 0, coinfected))

  # ── Stacked squares: one square per sample, youngest at the bottom ───────
  stack_data <- combination_counts %>%
    select(id, Age_Category, n) %>%
    arrange(id, Age_Category) %>%
    uncount(n) %>%
    group_by(id) %>%
    mutate(y = row_number()) %>%
    ungroup() %>%
    mutate(x_num = as.numeric(id))

  max_stack  <- max(stack_data$y)
  n_age <- length(age_levels)
  panel_rows <- max(max_stack, n_age + 1)   # +1 row for the legend header

  p_bars <- ggplot(stack_data, aes(x = x_num, y = y, fill = Age_Category)) +
    geom_tile(width = 0.7, height = 0.92, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = AGE_COLOURS, drop = FALSE) +
    scale_x_continuous(limits = c(0.5, n_combos + 0.5), expand = c(0.02, 0.02)) +
    scale_y_continuous(limits = c(0.5, panel_rows + 0.5), expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "none",  plot.margin = margin(10, 15, 1, 10))

  # ── Legend panel: colour key + N/total + % ───────────────────────────────
  legend_df <- age_stats %>%
    filter(Age_Category %in% age_levels) %>%
    mutate(Age_Category = factor(Age_Category, levels = age_levels)) %>%
    arrange(Age_Category) %>%
    mutate(y = panel_rows - row_number(),
           ratio = paste0(coinfected, "/", total_samples),
           pct = sprintf("%.1f", 100 * coinfected / total_samples))

  y_lo <- min(legend_df$y) - 0.5
  y_hi <- max(legend_df$y) + 0.5

  p_legend <- ggplot(legend_df) +
    annotate("rect", xmin = 3.05, xmax = 4.85, ymin = y_lo, ymax = y_hi,
             fill = "#f2f2f2") +
    annotate("segment", x = 3.05, xend = 3.05, y = y_lo, yend = y_hi,
             colour = "grey60", linewidth = 0.5) +
    geom_tile(aes(x = 0.9, y = y, fill = Age_Category),
              width = 0.3, height = 0.8) +
    geom_text(aes(x = 2.95, y = y, label = Age_Category),
              hjust = 1, size = 4, fontface = "bold") +
    geom_text(aes(x = 3.6, y = y, label = ratio), size = 3.5, fontface = "bold") +
    geom_text(aes(x = 4.4, y = y, label = pct),   size = 3.5, fontface = "bold") +
    annotate("text", x = c(3.6, 4.4), y = panel_rows,
             label = c("N/total", "%"), size = 4, fontface = "bold") +
    scale_fill_manual(values = AGE_COLOURS) +
    scale_x_continuous(limits = c(0.5, 4.9), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, panel_rows + 0.5), expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(10, 5, 1, 0))

  # ── Dot plot (unchanged) ─────────────────────────────────────────────────
  n_pathogens <- length(pathogen_order)
  reversed_order <- rev(pathogen_order)
  pathogen_y <- setNames(seq_along(reversed_order), reversed_order)

  df_long_plot <- df_long_plot %>%
    mutate(y_num = pathogen_y[as.character(Pathogen)])

  segment_data <- df_long_plot %>%
    filter(Present == 1) %>%
    group_by(id) %>%
    summarise(y_min = min(y_num), y_max = max(y_num), .groups = "drop") %>%
    filter(y_max > y_min) %>%
    mutate(x_num = as.numeric(id))

  p_dots <- ggplot(df_long_plot, aes(x = as.numeric(id), y = y_num)) +
    geom_point(data = filter(df_long_plot, Present == 0), color = "#D3D3D3", size = 4) +
    geom_point(data = filter(df_long_plot, Present == 1), color = "#2E4057", size = 4) +
    geom_segment(data = segment_data, aes(x = x_num, xend = x_num, y = y_min, yend = y_max), color = "#2E4057", linewidth = 0.8, inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq_len(n_combos), labels = NULL, expand = c(0.02, 0.02), limits = c(0.5, n_combos + 0.5)) +
    scale_y_continuous(breaks = seq_along(reversed_order), labels = reversed_order, limits = c(0.5, n_pathogens + 0.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x  = element_blank(),
          axis.text.y  = element_text(size = 12),
          panel.grid.major = element_line(color = "#F0F0F0"),
          panel.grid.minor = element_blank(),
          legend.position  = "none",
          plot.margin  = margin(2, 15, 10, 10))

  # ── Combine: bars + legend on top, dots below (legend column fixed at 3 in)
  final_plot <- p_bars + p_legend + p_dots + plot_spacer() +
    plot_layout(ncol    = 2, widths  = unit(c(1, 3), c("null", "in")), heights = c(panel_rows, n_pathogens)) +
    plot_annotation(
      title = title, theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )

  attr(final_plot, "recommended_height") <- max(4.5, 1 + (panel_rows + n_pathogens) * 0.18)
  return(final_plot)
}

# ── Run ───────────────────────────────────────────────────────────────────────
output_dir <- "results/Figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (name in names(datasets)) {

  if (!name %in% TARGET_SHEETS) next

  df <- datasets[[name]]
  colnames(df)[1:4] <- c("sample", "Pathogen", "Source", "Age")

  df_long <- df %>%
    mutate(Age_Category = sapply(Age, categorize_age)) %>%
    distinct(sample, Pathogen, .keep_all = TRUE)

  plot_title <- switch(name,
                       "Sheet2" = "Muco-kids",
                       "Sheet3" = "SUDI",
                       "Sheet4" = "Surveillance",
                       name)

  double_plot <- create_double_plot(df_long, plot_title, dataset_name = name)

plot_widths <- c("Sheet2" = 14, "Sheet3" = 10, "Sheet4" = 20)
plot_heights <- c("Sheet2" = 5,  "Sheet3" = 5, "Sheet4" = 5)   # optional

  if (!is.null(double_plot)) {
    plot_width <- plot_widths[[name]]
    rec_height <- if (name %in% names(plot_heights)) plot_heights[[name]] else
                    attr(double_plot, "recommended_height")

    base <- file.path(output_dir, paste0("FigureNew1_", plot_title))

    ggsave(paste0(base, ".png"), plot = double_plot, width = plot_width, height = rec_height, dpi = 300, bg = "white")
    ggsave(paste0(base, ".pdf"), plot = double_plot, width = plot_width, height = rec_height, device = cairo_pdf, bg = "white")
    ggsave(paste0(base, ".svg"), plot = double_plot, width = plot_width, height = rec_height, device = svglite::svglite, bg = "white")
  }
}