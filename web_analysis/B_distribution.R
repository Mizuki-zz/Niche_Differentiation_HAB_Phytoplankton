# Load necessary packages
library(readxl)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(extrafont)
library(tibble)  # Add this line
# Load Arial font (ensure the system has it installed)
loadfonts(device = "win")

# Custom function: load and process data for a single season (exclude diagonal)
load_season_data <- function(season) {
  file_name <- paste0("InterMatrix_", season, "_results.xlsx")
  
  # Read data and convert to matrix
  beta_matrix <- read_excel(file_name) %>%
    column_to_rownames(var = "Species") %>%
    as.matrix()
  
  # Mild outlier handling (Winsorizing)
  upper_limit <- quantile(abs(beta_matrix[beta_matrix != 0]), 0.99, na.rm = TRUE)
  beta_matrix[abs(beta_matrix) > upper_limit] <- sign(beta_matrix[abs(beta_matrix) > upper_limit]) * upper_limit
  
  # Set diagonal to 0 (intraspecific interactions are not considered)
  diag(beta_matrix) <- 0
  
  # Convert to long format and exclude zero values
  data.frame(
    Season = season,
    Row = rep(rownames(beta_matrix), each = ncol(beta_matrix)),
    Col = rep(colnames(beta_matrix), times = nrow(beta_matrix)),
    Value = as.vector(beta_matrix),
    stringsAsFactors = FALSE
  ) %>%
    filter(Value != 0) %>%  # Exclude zero values
    filter(Row != Col)      # Exclude diagonal
}

# Load data for all seasons
seasons <- c("spring", "summer", "autumn")
all_data <- map_df(seasons, load_season_data) %>%
  filter(is.finite(Value))

# Set season colors
season_colors <- c(
  spring = "#ff6666",
  summer = "#6699cc", 
  autumn = "#cc9900"
)

# Create visualization of interaction coefficient distribution
p1 <- ggplot(all_data, aes(x = Value, fill = Season)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = season_colors) +
  labs(title = "Distribution of Interaction Coefficients",
       x = "Interaction Coefficient (β)",
       y = "Density") +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = "top",
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  facet_wrap(~Season, ncol = 1, scales = "free_y")

# Create visualization of absolute interaction coefficient distribution
p2 <- ggplot(all_data, aes(x = abs(Value), fill = Season)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = season_colors) +
  labs(title = "Distribution of Absolute Interaction Coefficients",
       x = "Absolute Interaction Coefficient (|β|)",
       y = "Density") +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = "top",
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  facet_wrap(~Season, ncol = 1, scales = "free_y")

# Create boxplot comparing three seasons
p3 <- ggplot(all_data, aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = season_colors) +
  labs(title = "Comparison of Interaction Coefficients Across Seasons",
       x = "Season",
       y = "Interaction Coefficient (β)") +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")

# Create boxplot for absolute values
p4 <- ggplot(all_data, aes(x = Season, y = abs(Value), fill = Season)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = season_colors) +
  labs(title = "Comparison of Absolute Interaction Coefficients Across Seasons",
       x = "Season", 
       y = "Absolute Interaction Coefficient (|β|)") +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# Calculate and display quantile statistics
quantile_stats <- all_data %>%
  group_by(Season) %>%
  summarise(
    N = n(),
    Mean = mean(Value),
    SD = sd(Value),
    Median = median(Value),
    Q25 = quantile(Value, 0.25),
    Q75 = quantile(Value, 0.75),
    Min = min(Value),
    Max = max(Value)
  )

print("Interaction Coefficient Statistics:")
print(quantile_stats)

# Save visualization results
ggsave("Interaction_Coefficients_Distribution.pdf", p1, 
       width = 8, height = 10, device = cairo_pdf)

ggsave("Absolute_Interaction_Coefficients_Distribution.pdf", p2, 
       width = 8, height = 10, device = cairo_pdf)

ggsave("Interaction_Coefficients_Boxplot.pdf", p3, 
       width = 8, height = 6, device = cairo_pdf)

ggsave("Absolute_Interaction_Coefficients_Boxplot.pdf", p4, 
       width = 8, height = 6, device = cairo_pdf)

# The following section is not executed due to errors and lack of necessity
# Create combined report
library(patchwork)
combined_plot <- (p3 + p4) / (p1 + p2) + 
  plot_annotation(title = "Seasonal Comparison of Phytoplankton Interaction Coefficients",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", family = "Arial")))

ggsave("Combined_Interaction_Analysis.pdf", combined_plot, 
       width = 16, height = 12, device = cairo_pdf)

# Output statistical summary
cat("\n=== INTERACTION COEFFICIENT SUMMARY ===\n")
for(season in seasons) {
  season_data <- all_data %>% filter(Season == season)
  cat(sprintf("\n%s (n = %d):\n", season, nrow(season_data)))
  cat(sprintf("  Range: [%.3f, %.3f]\n", min(season_data$Value), max(season_data$Value)))
  cat(sprintf("  Mean ± SD: %.3f ± %.3f\n", mean(season_data$Value), sd(season_data$Value)))
  cat(sprintf("  Median: %.3f\n", median(season_data$Value)))
  cat(sprintf("  Positive interactions: %d (%.1f%%)\n", 
              sum(season_data$Value > 0), 100 * mean(season_data$Value > 0)))
  cat(sprintf("  Negative interactions: %d (%.1f%%)\n",
              sum(season_data$Value < 0), 100 * mean(season_data$Value < 0)))
}