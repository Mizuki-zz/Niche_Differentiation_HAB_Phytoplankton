library(readxl)
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(scales)
library(tibble)  # Added to fix "not find function column_to_rownames"
library(purrr)

# Custom parameter settings
THRESHOLD <- 0  # Interaction coefficient threshold
NODE_SIZE <- 5     # Node size
LINE_WIDTHS <- c(0.5, 1.0, 1.5)  # Three levels of edge width
POS_COLORS <- c("#FFC107", "#FF9800", "#FF5722")  # Promotional colors
NEG_COLORS <- c("#4FC3F7", "#2196F3", "#0D47A1")  # Inhibitory colors

# Season-specific node fill colors (new addition)
SEASON_COLORS <- c(
  spring = "#ff6666",  # Spring: light red
  summer = "#6699cc",  # Summer: light blue
  autumn = "#cc9900"   # Autumn: golden yellow
)

# Read species name mapping table
species_mapping <- read_excel("species_Eng.xlsx")
chinese_to_latin <- setNames(species_mapping$Latin_name, species_mapping$Chinese_name)

# 1. Read all seasonal data and combine
seasons <- c("spring", "summer", "autumn")
all_weights <- numeric(0)

for(season in seasons){
  file_name <- paste0("InterMatrix_", season, "_results.xlsx")
  beta_matrix <- read_excel(file_name) %>%
    column_to_rownames(var = "Species") %>%
    as.matrix()
  
  # Convert Chinese names to Latin names
  rownames(beta_matrix) <- chinese_to_latin[rownames(beta_matrix)]
  colnames(beta_matrix) <- chinese_to_latin[colnames(beta_matrix)]
  
  diag(beta_matrix) <- 0
  beta_matrix[abs(beta_matrix) <= THRESHOLD] <- 0
  
  season_weights <- abs(beta_matrix[beta_matrix != 0])
  all_weights <- c(all_weights, season_weights)
}

# 2. Unified classification criteria (Winsorization of extreme values) – not needed for topology analysis
upper_bound <- quantile(all_weights, probs = 0.99, na.rm = TRUE)
trimmed_weights <- pmin(all_weights, upper_bound)

global_quantiles <- quantile(
  trimmed_weights, 
  probs = c(0.33, 0.67), 
  na.rm = TRUE
)

# 3. Function to process each season (add node color parameter)
process_season <- function(season, quantiles, node_fill) {  # Added node_fill parameter
  file_name <- paste0("InterMatrix_", season, "_results.xlsx")
  beta_matrix <- read_excel(file_name) %>%
    column_to_rownames(var = "Species") %>%
    as.matrix()
  
  # Convert Chinese names to Latin names
  rownames(beta_matrix) <- chinese_to_latin[rownames(beta_matrix)]
  colnames(beta_matrix) <- chinese_to_latin[colnames(beta_matrix)]
  
  diag(beta_matrix) <- 0
  filtered_matrix <- beta_matrix
  filtered_matrix[abs(filtered_matrix) <= THRESHOLD] <- 0
  
  valid_species <- colnames(filtered_matrix)[colSums(abs(filtered_matrix)) > 0]
  filtered_matrix <- filtered_matrix[valid_species, valid_species]
  
  edges <- data.frame()
  for (i in rownames(filtered_matrix)) {
    for (j in colnames(filtered_matrix)) {
      beta <- filtered_matrix[i, j]
      if (beta != 0) {
        edges <- rbind(edges, data.frame(
          from = j, 
          to = i,
          weight = beta,
          abs_weight = abs(beta)
        ))
      }
    }
  }
  
  edges$level <- case_when(
    edges$abs_weight <= quantiles[1] ~ "low",
    edges$abs_weight <= quantiles[2] ~ "medium",
    TRUE ~ "high"
  )
  
  edges$color <- ifelse(
    edges$weight > 0,
    case_when(
      edges$level == "low" ~ POS_COLORS[1],
      edges$level == "medium" ~ POS_COLORS[2],
      edges$level == "high" ~ POS_COLORS[3]
    ),
    case_when(
      edges$level == "low" ~ NEG_COLORS[1],
      edges$level == "medium" ~ NEG_COLORS[2],
      edges$level == "high" ~ NEG_COLORS[3]
    )
  )
  
  edges$width <- case_when(
    edges$level == "low" ~ LINE_WIDTHS[1],
    edges$level == "medium" ~ LINE_WIDTHS[2],
    edges$level == "high" ~ LINE_WIDTHS[3]
  )
  
  graph <- graph_from_data_frame(edges, directed = TRUE)
  E(graph)$color <- edges$color
  E(graph)$width <- edges$width
  
  set.seed(123)
  # Calculate node layout
  layout_df <- create_layout(graph, layout = "linear", circular = TRUE)
  
  # Calculate the angle (in radians) for each node
  layout_df$angle <- atan2(layout_df$y, layout_df$x) * 180 / pi
  
  # Determine label position and adjustment parameters: Latin name position
  label_distance <- 1.05  # Multiplier for label distance from the center (relative to node position)
  
  # Improved calculation of label angle and alignment
  # Calculate label angle – align with the tangent direction of the circle
  layout_df$label_angle <- layout_df$angle
  
  # For points on the left half (angles between 90 and 270 degrees), adjust the angle
  left_side <- layout_df$angle > 90 & layout_df$angle < 270
  layout_df$label_angle[left_side] <- layout_df$angle[left_side] + 180
  
  # Ensure angles are within 0-360 degrees
  layout_df$label_angle <- layout_df$label_angle %% 360
  
  # Set horizontal alignment
  layout_df$hjust <- ifelse(layout_df$angle > 90 & layout_df$angle < 270, 1, 0)
  
  # Set vertical alignment
  layout_df$vjust <- ifelse(layout_df$angle > 90 & layout_df$angle < 270, 0.5, 0.5)
  
  # Adjust positions of labels on the left half so that the end of the word is near the node
  layout_df$label_x <- layout_df$x * label_distance
  layout_df$label_y <- layout_df$y * label_distance
  
  p_circle <- ggraph(layout_df) +
    geom_edge_arc(
      aes(edge_color = I(color), edge_width = I(width)),
      strength = 0.2,
      arrow = arrow(length = unit(3, 'mm'), type = "closed"), # Adjusted arrow size
      end_cap = circle(3, 'mm'),
      start_cap = circle(3, 'mm'),
      alpha = 0.8
    ) +
    # Modify node fill using the season‑specific color passed as argument
    geom_node_point(
      size = NODE_SIZE, 
      color = "gray30", 
      fill = node_fill,  # Use season‑specific color
      shape = 21
    ) +
    # Add Latin node labels
    geom_node_text(
      aes(x = label_x, y = label_y, label = name, 
          angle = label_angle, hjust = hjust, vjust = vjust),
      size = 5,  # Adjust font size if needed
      family = "Arial",
      fontface = "italic"
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      text = element_text(family = "Arial"),  # Changed to Arial font
      aspect.ratio = 1
    ) +
    expand_limits(x = c(-2, 2), y = c(-2, 2))  # Expand canvas to accommodate labels
  
  output_file <- paste0("phytoplankton_network_", season, ".pdf")
  ggsave(output_file, p_circle, width = 16, height = 16, device = cairo_pdf)
  
  return(edges)
}

# 4. Process all seasons and pass season‑specific colors
all_season_edges <- map2(
  seasons, 
  SEASON_COLORS[seasons],  # Pass corresponding season colors
  ~process_season(.x, global_quantiles, .y)
)

# 5. Create unified legend (English version)
legend_data <- data.frame(
  Strength = rep(c("Strong", "Medium", "Weak"), times = 2),
  Type = rep(c("Promotion", "Inhibition"), each = 3),
  Color = c(POS_COLORS[3:1], NEG_COLORS[3:1]),
  Width = rep(rev(LINE_WIDTHS), times = 2)
)

legend_data$Strength <- factor(legend_data$Strength, levels = c("Strong", "Medium", "Weak"))
legend_data$Type <- factor(legend_data$Type, levels = c("Promotion", "Inhibition"))
legend_data$start_x <- 1
legend_data$end_x <- 3

legend_plot <- ggplot(legend_data, aes(x = start_x, y = Strength)) +
  geom_segment(
    aes(xend = end_x, yend = Strength, color = Color, linewidth = Width),
    arrow = arrow(length = unit(3, 'mm'), type = "closed"), # Adjusted arrow size
    lineend = "round"
  ) +
  facet_grid(Strength ~ Type, scales = "free", space = "free") +
  scale_color_identity() +
  scale_linewidth_identity() +
  scale_x_continuous(breaks = NULL) +
  labs(y = "Interaction Strength") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 12),  # Changed to Arial font
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 10),
    strip.text.y.left = element_text(angle = 0),
    panel.spacing = unit(0.5, "cm")
  )

ggsave("network_legend.pdf", legend_plot, width = 6, height = 6, device = cairo_pdf)

# 6. Create season node color legend (English version)
season_legend <- ggplot() +
  geom_point(
    data = data.frame(
      Season = c("Spring", "Summer", "Autumn"),
      Color = SEASON_COLORS
    ),
    aes(x = Season, y = 1, fill = Color),
    shape = 21, size = 8
  ) +
  scale_fill_identity() +
  labs(fill = "Season") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, family = "Arial"),  # Changed to Arial font
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  coord_fixed(ratio = 0.1)

ggsave("season_legend.pdf", season_legend, width = 6, height = 2, device = cairo_pdf)

# Output global quantile information
cat("Global classification criteria (after Winsorization):\n")
cat(sprintf("Weak: ≤ %.3f\nMedium: %.3f - %.3f\nStrong: > %.3f",
            global_quantiles[1],
            global_quantiles[1],
            global_quantiles[2],
            global_quantiles[2]))