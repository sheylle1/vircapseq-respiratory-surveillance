
library(tidyverse)
library(ggplot2)

data <- read.csv("data/processed/data_reads.csv", stringsAsFactors = FALSE)

# Clean and prepare the data
analysis_ready <- data %>%
  mutate(
    # Clean totalReads by removing commas and converting to numeric
    totalReads = as.numeric(gsub(",", "", totalReads)),
    
    # Create Strategy as factor with the order you specified
    # Strategy = factor(Strategy, 
    #                   levels = c("24_V3", "60_P1", "96_P1", "96_P2", "48_P2"),
    #                   labels = c("24 Samples\non V3", "60 Samples\non P1", 
    #                              "96 Samples\non P1", "96 Samples\non P2", "48 Samples\non P2"))

   Strategy = factor(Strategy, 
                      levels = c("48_P2", "96_P2", "60_P1", "96_P1", "24_V3"),
                      labels = c(  "48 Samples\non P2", "96 Samples\non P2", "60 Samples\non P1", "96 Samples\non P1", "24 Samples\non V3"))                              
  )

median_data <- analysis_ready %>%
  group_by(Strategy) %>%
  summarise(median_val = median(totalReads/1e6, na.rm = TRUE))

y_max <- max(analysis_ready$totalReads/1e6, na.rm = TRUE)

p_boxplot <- analysis_ready %>%
  ggplot(aes(x = Strategy, y = totalReads/1e6)) +
  geom_boxplot(fill = "lightyellow", width = 0.5) +
  geom_jitter(alpha = 0.3, width = 0.2) +
  # Add median labels at the top
  geom_text(data = median_data,
            aes(x = Strategy, y = y_max * 1.3,  # Position slightly above the max
                label = sprintf("Median: %.2fM", median_val)),
            vjust = 0, size = 3.5) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")
  ) +
  labs(y = "Total Reads (Millions)",
       x = "Sequencing Strategy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 15)),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")


# Save the boxplot
ggsave("figures/Supfig1.png", plot = p_boxplot, width = 10, height = 6,  dpi = 300)
ggsave("figures/Supfig1.pdf", plot = p_boxplot, width = 10, height = 6)
