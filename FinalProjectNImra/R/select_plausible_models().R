#' Plausible Model Selection with AIC and Stability Filtering
#'
#' Selects plausible and stable models using AIC and average stability thresholds.
#' Optionally removes near-duplicate models using Jaccard similarity.
#'
#' @param full_data_models List of models from multi_path_forward().
#' @param stability_df Tibble of stability scores from stability_resampling().
#' @param delta_AIC Numeric; maximum AIC difference from best model.
#' @param tau Numeric; minimum average stability threshold.
#' @param jaccard_threshold Numeric; models with Jaccard similarity >= this will be treated as duplicates (default 0.9).
#' @return List of plausible and stable models.
#' @export
select_plausible_models <- function(full_data_models, stability_df,
                                    delta_AIC = 2, tau = 0.6,
                                    jaccard_threshold = 0.9) {
  # Compute AIC filter
  aic_values <- sapply(full_data_models, function(m) m$AIC)
  min_aic <- min(aic_values)
  plausible_models <- full_data_models[aic_values <= (min_aic + delta_AIC)]

  # Compute mean stability for each model
  stability_vec <- setNames(stability_df$Stability, stability_df$Variable)
  for (i in seq_along(plausible_models)) {
    vars <- plausible_models[[i]]$variables
    avg_stab <- mean(stability_vec[vars], na.rm = TRUE)
    plausible_models[[i]]$avg_stability <- avg_stab
  }

  # Filter by average stability threshold
  plausible_models <- Filter(function(m) m$avg_stability >= tau, plausible_models)

  # --- OPTIONAL: Remove near-duplicates via Jaccard similarity ---
  jaccard_similarity <- function(set1, set2) {
    intersect_len <- length(intersect(set1, set2))
    union_len <- length(union(set1, set2))
    return(intersect_len / union_len)
  }

  unique_models <- list()
  for (m in plausible_models) {
    if (length(unique_models) == 0) {
      unique_models <- list(m)
    } else {
      sims <- sapply(unique_models, function(u) jaccard_similarity(u$variables, m$variables))
      if (all(sims < jaccard_threshold)) {
        unique_models <- c(unique_models, list(m))
      }
    }
  }

  # Print summary
  if (length(unique_models) == 0) {
    cat("No plausible models found under current thresholds.\n")
  } else {
    cat("Plausible, stable models (after duplicate removal):\n")
    for (i in seq_along(unique_models)) {
      cat("Model", i, ":",
          paste(unique_models[[i]]$variables, collapse = ", "),
          "| AIC =", round(unique_models[[i]]$AIC, 3),
          "| Avg. Stability =", round(unique_models[[i]]$avg_stability, 2), "\n")
    }
  }

  return(unique_models)
}
