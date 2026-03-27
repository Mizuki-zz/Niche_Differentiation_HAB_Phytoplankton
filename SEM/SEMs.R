# ============================================================================
# Piecewise Structural Equation Modeling (SEM) Analysis
# Uses stepwise regression to remove non-significant paths for a testable model
# ============================================================================

# 1. Load required packages ----------------------------------------------------------
library(readxl)
library(piecewiseSEM)
library(nlme)
library(lme4)
library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)

# 2. Data reading and preprocessing function --------------------------------------------------
prepare_data <- function(file_path, sheet_name, env_vars) {
  # Read data
  data <- read_excel(file_path, sheet = sheet_name)
  
  # Select required variables
  selected_vars <- c("SampleID", "abundance", "diversity", "N", "F", env_vars)
  data <- data[, selected_vars]
  
  # Missing value handling: remove rows with missing values
  data_complete <- na.omit(data)
  cat(sprintf("\n%s season: original samples %d, complete samples %d (removed %d)\n", 
              sheet_name, nrow(data), nrow(data_complete), 
              nrow(data) - nrow(data_complete)))
  
  # Log transformation for skewed distributions
  data_complete$log_abundance <- log(data_complete$abundance + 1)
  
  # Standardize all variables used for modeling
  vars_to_scale <- c("diversity", "N", "F", env_vars, "log_abundance")
  data_scaled <- data_complete
  data_scaled[, vars_to_scale] <- scale(data_complete[, vars_to_scale])
  
  return(list(
    original = data_complete,
    scaled = data_scaled,
    n_obs = nrow(data_complete)
  ))
}

# 3. Build simplified initial SEM model (only theoretical paths) -----------------------------
build_initial_sem <- function(data, env_vars, season_name) {
  
  cat(sprintf("\n\n========== Building initial theoretical model for %s season ==========\n", season_name))
  
  # Core theoretical paths, excluding quadratic terms and all interactions
  # This ensures degrees of freedom for model testing
  
  # Path 1: Environmental factors → N
  formula_N <- as.formula(paste("N ~", paste(env_vars, collapse = " + ")))
  model_N <- lm(formula_N, data = data)
  
  # Path 2: Environmental factors + N → F (keep N→F due to high correlation)
  formula_F <- as.formula(paste("F ~", paste(c(env_vars, "N"), collapse = " + ")))
  model_F <- lm(formula_F, data = data)
  
  # Path 3: Environmental factors + N + F → diversity
  formula_div <- as.formula(paste("diversity ~", paste(c(env_vars, "N", "F"), collapse = " + ")))
  model_div <- lm(formula_div, data = data)
  
  # Path 4: Environmental factors + N + F + diversity → log_abundance
  formula_abun <- as.formula(paste("log_abundance ~", 
                                   paste(c(env_vars, "N", "F", "diversity"), collapse = " + ")))
  model_abun <- lm(formula_abun, data = data)
  
  # Build piecewise SEM
  psem_model <- psem(
    model_N,
    model_F,
    model_div,
    model_abun,
    data = data
  )
  
  return(psem_model)
}

# 4. Function to iteratively remove non-significant paths ---------------------------------------------
stepwise_simplify <- function(initial_model, data, season_name, alpha = 0.1) {
  
  cat(sprintf("\n--- Starting stepwise simplification for %s season (remove paths with p > %.2f) ---\n", season_name, alpha))
  
  current_models <- list(
    N = initial_model[[1]],
    F = initial_model[[2]],
    diversity = initial_model[[3]],
    log_abundance = initial_model[[4]]
  )
  
  iteration <- 1
  max_iterations <- 20
  
  while(iteration <= max_iterations) {
    cat(sprintf("\n========== Simplification iteration %d ==========\n", iteration))
    
    # Rebuild current psem model
    current_psem <- psem(
      current_models$N,
      current_models$F,
      current_models$diversity,
      current_models$log_abundance,
      data = data
    )
    
    # Get current model summary
    current_summary <- summary(current_psem, .progressBar = FALSE)
    coef_table <- current_summary$coefficients
    
    # Identify non‑critical paths with largest p‑value > alpha
    # Critical paths: among endogenous variables (N, F, diversity) are not removed
    critical_predictors <- c("N", "F", "diversity")
    
    # Mark removable paths
    coef_table$can_remove <- !(coef_table$Predictor %in% critical_predictors)
    coef_table$can_remove <- coef_table$can_remove & (coef_table$P.Value > alpha)
    
    removable <- coef_table[coef_table$can_remove, ]
    
    if(nrow(removable) == 0) {
      cat("\n✓ All retained paths have p ≤ %.2f, simplification complete\n", alpha)
      break
    }
    
    # Find the path with the largest p‑value
    max_p_idx <- which.max(removable$P.Value)
    path_to_remove <- removable[max_p_idx, ]
    
    response_var <- as.character(path_to_remove$Response)
    predictor_var <- as.character(path_to_remove$Predictor)
    
    cat(sprintf("\n→ Removing path: %s → %s (p = %.4f)\n", 
                predictor_var, response_var, path_to_remove$P.Value))
    
    # Update the corresponding model
    current_model <- current_models[[response_var]]
    current_formula <- formula(current_model)
    
    # Remove the predictor from the formula
    new_formula <- update(current_formula, paste(". ~ . -", predictor_var))
    
    # Refit the model
    current_models[[response_var]] <- lm(new_formula, data = data)
    
    cat(sprintf("  Updated formula: %s\n", deparse(new_formula)))
    
    iteration <- iteration + 1
  }
  
  if(iteration > max_iterations) {
    cat("\n⚠ Maximum iterations reached\n")
  }
  
  # Return the simplified model
  final_psem <- psem(
    current_models$N,
    current_models$F,
    current_models$diversity,
    current_models$log_abundance,
    data = data
  )
  
  return(final_psem)
}

# 5. Add significant missing paths (if needed) -----------------------------------------
add_missing_paths <- function(psem_model, data, season_name, alpha = 0.05) {
  
  cat(sprintf("\n--- Checking if missing paths need to be added for %s season ---\n", season_name))
  
  current_models <- list(
    N = psem_model[[1]],
    F = psem_model[[2]],
    diversity = psem_model[[3]],
    log_abundance = psem_model[[4]]
  )
  
  iteration <- 1
  max_iterations <- 5
  
  while(iteration <= max_iterations) {
    cat(sprintf("\n--- Path addition iteration %d ---\n", iteration))
    
    # Rebuild current model
    current_psem <- psem(
      current_models$N,
      current_models$F,
      current_models$diversity,
      current_models$log_abundance,
      data = data
    )
    
    # Check d‑separation
    dsep_tests <- dSep(current_psem, .progressBar = FALSE)
    
    if(is.null(dsep_tests) || nrow(dsep_tests) == 0) {
      cat("✓ No missing paths to test\n")
      break
    }
    
    # Identify significant missing paths
    significant_missing <- dsep_tests[dsep_tests$P.Value < alpha, ]
    
    if(nrow(significant_missing) == 0) {
      cat(sprintf("✓ All missing paths have p > %.2f, no need to add\n", alpha))
      break
    }
    
    # Select the path with the smallest p‑value
    min_p_idx <- which.min(significant_missing$P.Value)
    path_to_add <- significant_missing[min_p_idx, ]
    
    predictor <- as.character(path_to_add$Predictor)
    response <- as.character(path_to_add$Response)
    
    cat(sprintf("\n→ Adding path: %s → %s (p = %.4f)\n", 
                predictor, response, path_to_add$P.Value))
    
    # Update the model
    current_model <- current_models[[response]]
    current_formula <- formula(current_model)
    
    # Check if the predictor already exists
    if(!grepl(predictor, deparse(current_formula), fixed = TRUE)) {
      new_formula <- update(current_formula, paste(". ~ . +", predictor))
      current_models[[response]] <- lm(new_formula, data = data)
      cat(sprintf("  Updated formula: %s\n", deparse(new_formula)))
    } else {
      cat("  Path already exists, skipping\n")
      break
    }
    
    iteration <- iteration + 1
  }
  
  # Return the final model
  final_psem <- psem(
    current_models$N,
    current_models$F,
    current_models$diversity,
    current_models$log_abundance,
    data = data
  )
  
  return(final_psem)
}

# 6. Function to extract results ----------------------------------------------------------
extract_results <- function(psem_model, season_name) {
  
  # Get model summary
  model_summary <- summary(psem_model, .progressBar = FALSE)
  
  # 1. Path coefficients
  coef_table <- as.data.frame(model_summary$coefficients)
  coef_table$Season <- season_name
  coef_table$Path <- paste(coef_table$Predictor, "→", coef_table$Response)
  
  # Add significance markers
  coef_table$Significance <- ifelse(coef_table$P.Value < 0.001, "***",
                                    ifelse(coef_table$P.Value < 0.01, "**",
                                           ifelse(coef_table$P.Value < 0.05, "*", "ns")))
  
  # Reorder columns
  coef_table <- coef_table[, c("Season", "Path", "Response", "Predictor", 
                               "Estimate", "Std.Error", "DF", "Crit.Value", 
                               "P.Value", "Std.Estimate", "Significance")]
  
  # 2. Model fit indices
  fisher_c <- model_summary$Cstat$Fisher.C
  df <- model_summary$Cstat$df
  p_value <- model_summary$Cstat$P.Value
  
  # Handle NA values
  if(is.na(fisher_c)) fisher_c <- 0
  if(is.na(df)) df <- 0
  if(is.na(p_value)) p_value <- 1
  
  fit_stats <- data.frame(
    Season = season_name,
    Fisher_C = fisher_c,
    df = df,
    P_value = p_value,
    AIC = model_summary$AIC,
    stringsAsFactors = FALSE
  )
  
  # Add fit assessment
  if(df > 0) {
    fit_stats$Model_Fit <- ifelse(p_value > 0.05, "Good", "Poor")
    fit_stats$Chi_sq_test <- sprintf("χ²(%.0f) = %.2f, p = %.4f", df, fisher_c, p_value)
  } else {
    fit_stats$Model_Fit <- "Saturated (df=0)"
    fit_stats$Chi_sq_test <- "Not testable"
  }
  
  # 3. R² values
  r2_list <- list()
  for(i in seq_along(psem_model)) {
    model <- psem_model[[i]]
    if(inherits(model, "lm")) {
      response_var <- all.vars(formula(model))[1]
      r2_val <- summary(model)$r.squared
      adj_r2_val <- summary(model)$adj.r.squared
      
      r2_list[[length(r2_list) + 1]] <- data.frame(
        Season = season_name,
        Variable = response_var,
        R_squared = round(r2_val, 4),
        Adj_R_squared = round(adj_r2_val, 4),
        stringsAsFactors = FALSE
      )
    }
  }
  r2_table <- do.call(rbind, r2_list)
  
  return(list(
    coefficients = coef_table,
    fit_indices = fit_stats,
    r_squared = r2_table,
    summary = model_summary
  ))
}

# 7. Visualization function ------------------------------------------------------------
plot_sem_results <- function(results_list, output_dir = ".") {
  
  # Plot comparison of significant paths
  all_coefs <- do.call(rbind, lapply(results_list, function(x) x$coefficients))
  sig_coefs <- all_coefs[all_coefs$P.Value < 0.05, ]
  
  if(nrow(sig_coefs) > 0) {
    p1 <- ggplot(sig_coefs, aes(x = reorder(Path, Std.Estimate), 
                                y = Std.Estimate, fill = Season)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      geom_errorbar(aes(ymin = Std.Estimate - Std.Error, 
                        ymax = Std.Estimate + Std.Error),
                    position = position_dodge(0.7), width = 0.3) +
      coord_flip() +
      theme_minimal() +
      labs(title = "Comparison of Standardized Coefficients for Significant Paths (p < 0.05)",
           x = "Path", y = "Standardized Coefficient") +
      theme(axis.text.y = element_text(size = 7),
            legend.position = "bottom")
    
    ggsave(file.path(output_dir, "significant_paths_comparison.png"), 
           p1, width = 14, height = 10, dpi = 300)
    cat("\n✓ Path coefficient comparison plot saved\n")
  }
  
  # R² value comparison
  all_r2 <- do.call(rbind, lapply(results_list, function(x) x$r_squared))
  
  p2 <- ggplot(all_r2, aes(x = Variable, y = R_squared, fill = Season)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = sprintf("%.3f", R_squared)), 
              position = position_dodge(0.9), vjust = -0.5, size = 3) +
    theme_minimal() +
    labs(title = "Comparison of Variance Explained (R²) by Variable",
         x = "Variable", y = "R²") +
    ylim(0, max(all_r2$R_squared) * 1.15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
  
  ggsave(file.path(output_dir, "r_squared_comparison.png"), 
         p2, width = 10, height = 6, dpi = 300)
  cat("✓ R² comparison plot saved\n")
  
  # Model fit comparison
  all_fits <- do.call(rbind, lapply(results_list, function(x) x$fit_indices))
  
  p3 <- ggplot(all_fits[all_fits$df > 0, ], 
               aes(x = Season, y = P_value)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    geom_text(aes(label = sprintf("p = %.3f", P_value)), 
              vjust = -0.5, size = 4) +
    theme_minimal() +
    labs(title = "Model Fit (Fisher's C Test)",
         subtitle = "p > 0.05 indicates good fit",
         x = "Season", y = "P-value") +
    ylim(0, max(c(all_fits$P_value[all_fits$df > 0], 0.1)) * 1.2)
  
  ggsave(file.path(output_dir, "model_fit_comparison.png"), 
         p3, width = 8, height = 6, dpi = 300)
  cat("✓ Model fit comparison plot saved\n")
}

# ============================================================================
# Main program
# ============================================================================

main <- function() {
  
  # Set file path
  file_path <- "SEMs_3.xlsx"
  
  # Define environmental variables for each season
  env_config <- list(
    spring = c("Temperature", "Salinity"),
    summer = c("Salinity", "Oxygen"),
    autumn = c("Oxygen", "pH")
  )
  
  # Store all results
  all_results <- list()
  all_coef_tables <- list()
  all_fit_tables <- list()
  all_r2_tables <- list()
  
  # Model each of the three seasons
  for(season in c("spring", "summer", "autumn")) {
    
    cat(sprintf("\n\n%s\n", paste(rep("=", 80), collapse = "")))
    cat(sprintf("Processing %s season data\n", toupper(season)))
    cat(sprintf("%s\n", paste(rep("=", 80), collapse = "")))
    
    # 1. Data preparation
    data_prep <- prepare_data(file_path, season, env_config[[season]])
    
    # 2. Build initial simplified model
    initial_model <- build_initial_sem(
      data_prep$scaled, 
      env_config[[season]], 
      season
    )
    
    cat("\n---------- Initial model summary ----------\n")
    initial_summary <- summary(initial_model, .progressBar = FALSE)
    print(initial_summary)
    
    # 3. Iteratively remove non‑significant paths
    simplified_model <- stepwise_simplify(
      initial_model, 
      data_prep$scaled, 
      season,
      alpha = 0.1  # Use a more lenient threshold
    )
    
    cat("\n---------- Simplified model summary ----------\n")
    simplified_summary <- summary(simplified_model, .progressBar = FALSE)
    print(simplified_summary)
    
    # 4. Add significant missing paths (if any)
    final_model <- add_missing_paths(
      simplified_model,
      data_prep$scaled,
      season,
      alpha = 0.05
    )
    
    # 5. Print final model
    cat(sprintf("\n\n%s\n", paste(rep("-", 80), collapse = "")))
    cat("Final optimized model summary\n")
    cat(sprintf("%s\n", paste(rep("-", 80), collapse = "")))
    final_summary <- summary(final_model, .progressBar = FALSE)
    print(final_summary)
    
    # 6. Extract results
    results <- extract_results(final_model, season)
    all_results[[season]] <- results
    
    all_coef_tables[[season]] <- results$coefficients
    all_fit_tables[[season]] <- results$fit_indices
    all_r2_tables[[season]] <- results$r_squared
    
    # 7. Print key metrics
    cat(sprintf("\n\n%s\n", paste(rep("=", 80), collapse = "")))
    cat(sprintf("Key results for %s season\n", toupper(season)))
    cat(sprintf("%s\n", paste(rep("=", 80), collapse = "")))
    
    cat(sprintf("Model fit indices:\n"))
    cat(sprintf("  Fisher's C: %.4f\n", results$fit_indices$Fisher_C))
    cat(sprintf("  Degrees of freedom (df): %d\n", results$fit_indices$df))
    cat(sprintf("  P-value: %.4f\n", results$fit_indices$P_value))
    cat(sprintf("  AIC: %.4f\n", results$fit_indices$AIC))
    cat(sprintf("  Model fit: %s\n", results$fit_indices$Model_Fit))
    
    cat(sprintf("\nR² values for each variable:\n"))
    print(results$r_squared, row.names = FALSE)
    
    abundance_r2 <- results$r_squared$R_squared[results$r_squared$Variable == "log_abundance"]
    cat(sprintf("\n★ R² for log_abundance: %.4f\n", abundance_r2))
  }
  
  # 8. Combine and export results
  cat(sprintf("\n\n%s\n", paste(rep("=", 80), collapse = "")))
  cat("Combine results and export to Excel file\n")
  cat(sprintf("%s\n", paste(rep("=", 80), collapse = "")))
  
  combined_coefs <- do.call(rbind, all_coef_tables)
  combined_fits <- do.call(rbind, all_fit_tables)
  combined_r2 <- do.call(rbind, all_r2_tables)
  
  output_file <- "SEM_Results_Summary.xlsx"
  
  write_xlsx(
    list(
      "Path_Coefficients" = combined_coefs,
      "Model_Fit_Indices" = combined_fits,
      "R_Squared" = combined_r2
    ),
    path = output_file
  )
  
  cat(sprintf("\n✓ Results saved to: %s\n", output_file))
  
  # 9. Visualization
  cat("\nGenerating visualizations...\n")
  plot_sem_results(all_results)
  
  # 10. Summary of significant paths
  cat(sprintf("\n\n%s\n", paste(rep("=", 80), collapse = "")))
  cat("Summary of significant paths (p < 0.05)\n")
  cat(sprintf("%s\n", paste(rep("=", 80), collapse = "")))
  
  for(season in c("spring", "summer", "autumn")) {
    cat(sprintf("\n【%s Season】\n", toupper(season)))
    coefs <- all_results[[season]]$coefficients
    sig_paths <- coefs[coefs$P.Value < 0.05, 
                       c("Path", "Std.Estimate", "P.Value", "Significance")]
    sig_paths <- sig_paths[order(sig_paths$P.Value), ]
    
    if(nrow(sig_paths) > 0) {
      sig_paths$Std.Estimate <- round(sig_paths$Std.Estimate, 4)
      sig_paths$P.Value <- sprintf("%.4f", sig_paths$P.Value)
      print(sig_paths, row.names = FALSE)
    } else {
      cat("  No significant paths\n")
    }
  }
  
  cat(sprintf("\n\n%s\n", paste(rep("=", 80), collapse = "")))
  cat("Analysis complete!\n")
  cat(sprintf("%s\n", paste(rep("=", 80), collapse = "")))
  
  return(all_results)
}

# Run main program
results <- main()

cat("\nProgram execution finished!\n")