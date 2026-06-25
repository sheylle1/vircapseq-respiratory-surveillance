library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readxl)
library(reshape2)

sheets <- excel_sheets("/data/processed/data_coinfect.xlsx")
datasets <- list()
for (i in 1:length(sheets)) {
  datasets[[i]] <- read_excel("/data/processed/data_coinfect.xlsx", sheet = sheets[i])
  names(datasets)[i] <- sheets[i]
}

categorize_age <- function(age) {
  if (is.na(age)) return("Unknown")
  age <- as.numeric(age)
  if (age < 1) return("<1 years")
  else if (age >= 1 & age < 5) return("1-5 years")
  else if (age >= 5 & age < 18) return("5-18 years")
  else if (age >= 18 & age <= 60) return("18-60 years")
  else if (age > 60) return(">60 years")
  else return("Unknown")
}

create_double_plot <- function(df_long, title, dataset_name = NULL) {
  binary_matrix <- df_long %>%
    mutate(Present = 1) %>%
    pivot_wider(id_cols = c(sample, Age_Category), 
                names_from = Pathogen, 
                values_from = Present, 
                values_fill = list(Present = 0)) %>%
    as.data.frame()
  
  pathogen_cols <- setdiff(colnames(binary_matrix), c("sample", "Age_Category"))
  
  if (length(pathogen_cols) == 0) {
    return(NULL)
  }
  
  binary_matrix$CombinationID <- apply(binary_matrix[, pathogen_cols, drop = FALSE], 1, 
                                       function(x) paste(x, collapse = ""))
  
  combo_totals <- binary_matrix %>%
    group_by(CombinationID) %>%
    summarise(total_n = n(), .groups = 'drop') %>%
    arrange(desc(total_n))
  
  combination_counts <- binary_matrix %>%
    group_by(CombinationID, Age_Category) %>%
    summarise(n = n(), .groups = 'drop') %>%
    left_join(
      binary_matrix %>%
        select(-sample, -Age_Category) %>%
        distinct(CombinationID, .keep_all = TRUE),
      by = "CombinationID"
    )
  
  combination_counts <- combination_counts %>%
    rowwise() %>%
    mutate(num_pathogens = sum(c_across(all_of(pathogen_cols)))) %>%
    ungroup() %>%
    filter(num_pathogens >= 2)
  
  if (nrow(combination_counts) == 0) {
    return(NULL)
  }
  
  combination_counts <- combination_counts %>%
    left_join(combo_totals, by = "CombinationID") %>%
    arrange(desc(total_n), CombinationID)
  
  combo_labels <- combination_counts %>%
    distinct(CombinationID, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(
      Pathogens = paste(pathogen_cols[c_across(all_of(pathogen_cols)) == 1], collapse = " + ")
    ) %>%
    ungroup() %>%
    arrange(desc(total_n), CombinationID) %>%
    mutate(id = factor(as.character(seq_len(n())), levels = as.character(seq_len(n())))) %>%
    select(id, CombinationID, Pathogens, total_n)
  
  combo_id_map <- combo_labels %>% select(CombinationID, id)
  
  combination_counts <- combination_counts %>%
    left_join(combo_id_map, by = "CombinationID")
  
  pathogens_present <- combination_counts %>%
    select(all_of(pathogen_cols)) %>%
    summarise(across(everything(), ~any(. == 1))) %>%
    pivot_longer(everything(), names_to = "Pathogen", values_to = "present") %>%
    filter(present == TRUE) %>%
    pull(Pathogen)
  
  df_long_plot <- combination_counts %>%
    pivot_longer(cols = all_of(pathogen_cols),
                 names_to = "Pathogen",
                 values_to = "Present") %>%
    filter(Pathogen %in% pathogens_present)
  
  pathogen_order <- sort(unique(pathogens_present))
  
  df_long_plot$Pathogen <- factor(df_long_plot$Pathogen, levels = pathogen_order)
  
  combo_order <- combo_labels$id
  combination_counts$id <- factor(combination_counts$id, levels = combo_order)
  df_long_plot$id <- factor(df_long_plot$id, levels = combo_order)
  
  age_levels <- c("<1 years", "1-5 years", "5-18 years", "18-60 years", ">60 years", "Unknown")
  age_levels <- age_levels[age_levels %in% unique(combination_counts$Age_Category)]
  
  heatmap_data <- expand.grid(
    id = factor(combo_order, levels = combo_order),
    Age_Category = factor(age_levels, levels = age_levels),
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      combination_counts %>% select(id, Age_Category, n),
      by = c("id", "Age_Category")
    ) %>%
    mutate(n = ifelse(is.na(n), 0, n))
  
  max_count <- max(heatmap_data$n, na.rm = TRUE)
  if (max_count == 0) max_count <- 1
  
  is_surveillance <- !is.null(dataset_name) && dataset_name == "Sheet4"
  
  if (is_surveillance) {
    heatmap_data$text_color <- ifelse(heatmap_data$n <= 10, "black", "white")
  } else {
    heatmap_data$text_color <- "white"
  }
  
  p_heatmap <- ggplot(heatmap_data, aes(x = id, y = Age_Category, fill = n)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(n > 0, as.character(n), ""), 
                  color = text_color), 
              size = 4, fontface = "bold", show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradientn(
      colors = c("#FFFFFF", "#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", 
                 "#D95F0E", "#993404"),
      name = "Count",
      limits = c(0, max_count),
      na.value = "white"
    ) +
    scale_y_discrete(limits = rev(age_levels)) +
    scale_x_discrete(labels = NULL, expand = c(0, 0)) +
    labs(x = NULL, y = NULL, fill = "Count") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.4, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.margin = margin(5, 15, 1, 10),
      plot.title.position = "plot"
    ) 
  
  p_dots <- ggplot(df_long_plot, aes(x = id, y = Pathogen, color = factor(Present))) +
    geom_point(size = 4) +
    scale_color_manual(values = c("0" = "#D3D3D3", "1" = "#2E4057"), guide = "none") +
    scale_y_discrete(limits = rev) +
    scale_x_discrete(labels = NULL) +
    labs(x = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_blank(),
      panel.grid.major = element_line(color = "#F0F0F0"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 15, 10, 10)
    ) 
  
  if (is_surveillance) {
    final_plot <- (p_heatmap / p_dots) + 
      plot_layout(heights = c(0.04, 0.1), guides = "collect") &
      theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10)
      ) &
      plot_annotation(
        title = paste(title),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
      )
    
    n_combos <- length(combo_order)
    rec_height <- max(3, 1 + n_combos * 0.1)
  } else {
    final_plot <- (p_heatmap / p_dots) + 
      plot_layout(heights = c(0.15, 0.8), guides = "collect") &
      theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10)
      ) &
      plot_annotation(
        title = paste(title),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
      )
    
    n_combos <- length(combo_order)
    rec_height <- max(3, 1 + n_combos * 0.1)
  }
  
  attr(final_plot, "recommended_height") <- rec_height
  
  return(final_plot)
}

for (name in names(datasets)) {
  df <- datasets[[name]]
  colnames(df)[1:4] <- c("sample", "Pathogen", "Source", "Age")
  df_long <- df
  
  df_long <- df_long %>%
    mutate(Age_Category = sapply(Age, categorize_age)) %>%
    distinct(sample, Pathogen, .keep_all = TRUE)
  
  plot_title <- switch(name,
                       "Sheet2" = "Mucokids",
                       "Sheet3" = "SUDI",
                       "Sheet4" = "Surveillance",
                       name)
  
  double_plot <- create_double_plot(df_long, plot_title, dataset_name = name)
  if (!is.null(double_plot)) {
    rec_height <- attr(double_plot, "recommended_height")
    if (is.null(rec_height)) rec_height <- 9
    ggsave(paste0("double_plot_", name, ".png"), double_plot, 
           width = 16, height = rec_height, dpi = 300)
  }
}