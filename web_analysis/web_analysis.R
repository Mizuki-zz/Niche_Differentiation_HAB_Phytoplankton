# Load necessary packages
library(readxl)
library(igraph)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(purrr)
library(tidyr)
library(tibble)

# Custom parameter settings
THRESHOLD <- 0  # Interaction coefficient threshold, change as needed.
N_SIMULATIONS <- 50 # Number of simulations for robustness analysis

# 1. Modified network topology analysis function
analyze_topology <- function(graph) {
  # Create absolute weight graph for path calculation
  abs_graph <- set_edge_attr(graph, "abs_weight", value = abs(E(graph)$weight))
  
  # Calculate basic topological metrics
  metrics <- list(
    Nodes = vcount(graph),
    Edges = ecount(graph),
    Density = edge_density(graph),
    Mean_degree = mean(degree(graph, mode = "all")),
    Transitivity = transitivity(graph, type = "global"),
    Diameter = diameter(abs_graph, directed = TRUE, weights = E(abs_graph)$abs_weight),
    Avg_path_length = mean_distance(abs_graph, directed = TRUE, weights = E(abs_graph)$abs_weight)
  )
  
  # Calculate modularity using absolute weights
  # Create undirected graph for community detection (walktrap algorithm requires undirected graph)
  undirected_abs_graph <- as.undirected(
    abs_graph, 
    mode = "each", 
    edge.attr.comb = list(weight = "sum", abs_weight = "sum")
  )
  
  # Perform community detection using absolute weights
  comm <- cluster_walktrap(undirected_abs_graph, weights = E(undirected_abs_graph)$abs_weight)
  metrics$Modularity <- modularity(comm)
  
  # Calculate node centrality metrics
  centrality <- list(
    Degree_centrality = degree(graph, mode = "all"),
    Betweenness_centrality = betweenness(abs_graph, directed = TRUE, weights = E(abs_graph)$abs_weight),
    Closeness_centrality = closeness(abs_graph, mode = "all", weights = E(abs_graph)$abs_weight),
    Eigenvector_centrality = eigen_centrality(abs_graph)$vector
  )
  
  return(list(topology_metrics = metrics, centrality = centrality))
}

# 2. Interaction strength and direction analysis function
analyze_interactions <- function(edges) {
  # Basic statistics
  interaction_stats <- list(
    Total_interactions = nrow(edges),
    Positive_interactions = sum(edges$weight > 0),
    Negative_interactions = sum(edges$weight < 0),
    Prop_positive = mean(edges$weight > 0),
    Mean_strength = mean(abs(edges$weight)),
    SD_strength = sd(abs(edges$weight))
  )
  
  # Strength distribution
  strength_dist <- density(abs(edges$weight))
  
  return(list(stats = interaction_stats, distribution = strength_dist))
}

# 3. Network stability analysis function
analyze_stability <- function(graph) {
  # Create absolute weight graph for stability calculation
  abs_graph <- set_edge_attr(graph, "abs_weight", value = abs(E(graph)$weight))
  
  # Natural connectivity using absolute weights
  natural_connectivity <- function(g) {
    adj_matrix <- as.matrix(get.adjacency(g, attr = "abs_weight"))
    eigenvalues <- eigen(adj_matrix, only.values = TRUE)$values
    lambda <- max(Re(eigenvalues))
    return(lambda)
  }
  
  # Simulate random node attacks (using unweighted graph)
  simulate_node_attack <- function(g) {
    nodes <- vcount(g)
    robustness <- numeric(nodes)
    
    # Create unweighted copy for attack simulation
    unweighted_g <- delete_edge_attr(g, "weight")
    if ("abs_weight" %in% edge_attr_names(g)) {
      unweighted_g <- delete_edge_attr(unweighted_g, "abs_weight")
    }
    
    for (i in 1:nodes) {
      nodes_to_remove <- sample(1:nodes, i)
      g_attacked <- delete_vertices(unweighted_g, nodes_to_remove)
      
      components <- components(g_attacked, mode = "weak")
      robustness[i] <- max(components$csize) / vcount(g_attacked)
    }
    return(robustness)
  }
  
  # Average over multiple simulations
  robustness_matrix <- replicate(N_SIMULATIONS, simulate_node_attack(graph))
  avg_robustness <- rowMeans(robustness_matrix)
  
  return(list(
    natural_connectivity = natural_connectivity(abs_graph),
    robustness_curve = avg_robustness
  ))
}

# 4. Main analysis function
analyze_season <- function(season) {
  # Read data
  file_name <- paste0("InterMatrix_", season, "_results.xlsx")
  beta_matrix <- read_excel(file_name) %>%
    column_to_rownames(var = "Species") %>%
    as.matrix()
  
  # Filter weak interactions
  diag(beta_matrix) <- 0
  beta_matrix[abs(beta_matrix) <= THRESHOLD] <- 0
  
  # Build network
  valid_species <- colnames(beta_matrix)[colSums(abs(beta_matrix)) > 0]
  filtered_matrix <- beta_matrix[valid_species, valid_species]
  
  # Create edge list
  edges <- data.frame()
  for (i in rownames(filtered_matrix)) {
    for (j in colnames(filtered_matrix)) {
      beta <- filtered_matrix[i, j]
      if (beta != 0) {
        edges <- rbind(edges, data.frame(
          from = j, 
          to = i,
          weight = beta,
          abs_weight = abs(beta),
          Season = season
        ))
      }
    }
  }
  
  # Create graph object
  graph <- graph_from_data_frame(edges, directed = TRUE)
  
  # Perform analysis
  topology_results <- analyze_topology(graph)
  interaction_results <- analyze_interactions(edges)
  stability_results <- analyze_stability(graph)
  
  # Prepare output data
  list(
    season = season,
    edges = edges,
    graph = graph,
    topology_metrics = topology_results$topology_metrics,
    centrality = topology_results$centrality,
    interaction_stats = interaction_results$stats,
    strength_distribution = interaction_results$distribution,
    stability = stability_results
  )
}

# 5. Analyze all seasons
seasons <- c("spring", "summer", "autumn")
all_results <- map(seasons, analyze_season)
names(all_results) <- seasons

# 6. Result summary and output
wb <- createWorkbook()

# Summarize topology metrics
topology_summary <- map_dfr(all_results, function(res) {
  as.data.frame(res$topology_metrics)
}, .id = "Season")

addWorksheet(wb, "Topology_Metrics")
writeData(wb, "Topology_Metrics", topology_summary)

# Summarize interaction statistics
interaction_summary <- map_dfr(all_results, function(res) {
  as.data.frame(res$interaction_stats)
}, .id = "Season")

addWorksheet(wb, "Interaction_Stats")
writeData(wb, "Interaction_Stats", interaction_summary)

# Function to calculate Composite_Score for all species
calculate_all_composite_scores <- function(centrality_list) {
  centrality_df <- data.frame(
    Species = names(centrality_list$Degree_centrality),
    Degree_centrality = centrality_list$Degree_centrality,
    Betweenness_centrality = centrality_list$Betweenness_centrality,
    Closeness_centrality = centrality_list$Closeness_centrality,
    Eigenvector_centrality = centrality_list$Eigenvector_centrality
  )
  
  # Normalize centrality metrics and calculate composite score
  centrality_df %>%
    mutate(across(where(is.numeric), scales::rescale)) %>%
    mutate(Composite_Score = rowMeans(select(., -Species))) %>%
    arrange(desc(Composite_Score))
}

# Centrality results (including Composite_Score for all species)
for (season in seasons) {
  centrality_df <- calculate_all_composite_scores(all_results[[season]]$centrality)
  
  sheet_name <- paste0("Centrality_", season)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, centrality_df)
}

# 7. Visualization
# 7.1 Topology metrics comparison
topology_plot <- topology_summary %>%
  select(Season, Nodes, Edges, Density, Mean_degree, Transitivity, Modularity) %>%
  pivot_longer(-Season, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Season, y = Value, fill = Season)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Network Topology Metrics Comparison",
       y = "Value") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("Topology_Metrics_Comparison.pdf", topology_plot, width = 10, height = 8)

# 7.2 Interaction strength distribution
strength_plot <- bind_rows(map(all_results, ~as.data.frame(.x$edges))) %>%
  ggplot(aes(x = abs_weight, fill = Season)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Season, ncol = 1) +
  labs(title = "Interaction Strength Distribution",
       x = "Absolute Interaction Strength",
       y = "Density") +
  theme_minimal()

ggsave("Interaction_Strength_Distribution.pdf", strength_plot, width = 8, height = 10)

# 7.3 Robustness curve
robustness_data <- map_dfr(all_results, function(res) {
  data.frame(
    Season = res$season,
    Nodes_removed = 1:length(res$stability$robustness_curve),
    Robustness = res$stability$robustness_curve
  )
})

robustness_plot <- robustness_data %>%
  ggplot(aes(x = Nodes_removed, y = Robustness, color = Season)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Network Robustness to Random Node Attacks",
       x = "Number of Nodes Removed",
       y = "Proportion of Largest Connected Component") +
  theme_minimal()

ggsave("Network_Robustness.pdf", robustness_plot, width = 8, height = 6)

# 8. Key species identification (top 5 based on Composite_Score)
identify_key_species <- function(centrality_df_with_composite) {
  centrality_df_with_composite %>%
    arrange(desc(Composite_Score)) %>%
    head(5)  # Take top 5 key species per season
}

key_species_summary <- map_dfr(seasons, function(season) {
  centrality_df <- calculate_all_composite_scores(all_results[[season]]$centrality)
  key_species <- identify_key_species(centrality_df) %>%
    mutate(Season = season, .before = Species)
  return(key_species)
})

# Save key species results
addWorksheet(wb, "Key_Species")
writeData(wb, "Key_Species", key_species_summary)

# Save Excel file
saveWorkbook(wb, "Phytoplankton_Network_Analysis.xlsx", overwrite = TRUE)

# 9. Output natural connectivity results
cat("Natural Connectivity Results:\n")
for (season in seasons) {
  cat(sprintf("%s: %.4f\n", season, all_results[[season]]$stability$natural_connectivity))
}

# 10. Output analysis summary
cat("\nAnalysis Summary:\n")
cat(sprintf("- Analyzed %d seasons: %s\n", length(seasons), paste(seasons, collapse = ", ")))
cat(sprintf("- Threshold for significant interactions: %.2f\n", THRESHOLD))
cat("- Generated Excel file: Phytoplankton_Network_Analysis.xlsx\n")
cat("- Generated PDF visualizations: Topology_Metrics_Comparison.pdf, Interaction_Strength_Distribution.pdf, Network_Robustness.pdf\n")