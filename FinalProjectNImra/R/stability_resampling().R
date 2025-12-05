#' Stability Estimation with Bootstrap Resampling
#'
#' Estimates the stability of variable selection using repeated bootstrap samples.
#'
#' @param data Data frame.
#' @param response Response variable.
#' @param B Number of bootstrap samples.
#' @param K,epsilon,delta,L Parameters passed to multi_path_forward().
#' @param model_type Model type ("linear" or "logistic"). Automatically detected if not provided.
#' @return A tibble with variable names and their stability scores.
#' @export
stability <- function(X, response, B = 100, K = 5, epsilon = 1e-6,
                      delta = 2, L = 25, model_type = NULL) {

  result_matrix <- matrix(0, nrow = B, ncol = length(setdiff(names(X), response)),
                          dimnames = list(NULL, setdiff(names(X), response)))

  resampling <- vector("list", B)
  predictors <- colnames(result_matrix)
  p_length <- length(predictors)

  for (b in 1:B) {
    # Bootstrap
    index_b <- sample(1:nrow(X), size = nrow(X), replace = TRUE)
    X_b <- X[index_b, , drop = FALSE]

    # Fit model
    models_b <- build_paths(X_b, response,
                            K = K, epsilon = epsilon,
                            delta = delta, L = L,
                            model_type = model_type)

    # Feature inclusion: counts for # times each predictor was used
    count_predictor_usage <- function(models_b) {

      # extract all "vars" from all models across all frontiers
      all_vars <- unlist(
        lapply(models_b$path_forest$frontiers, function(level)
          lapply(level, function(model) model$vars)),
        recursive = TRUE
      )

      # tabulate
      usage <- table(all_vars)

      return(usage)
    }

    counts <- count_predictor_usage(models_b)
    counts

    # Proportions for Bootstrap sample
    z_j_b <- counts / (ifelse(length(models_b) > 0, length(models_b), 1)) # calculates ( z_j^{(b)} )
    result_matrix[b, ] <- z_j_b

    # Store proportions: is this redundant? isn't the line above this storing the bootstrap?
    resampling[[b]] <- z_j_b
  }

  # Calculate stability
  pi <- colMeans(result_matrix)
  path_stability <- data.frame(variable = predictors, pi = pi)
  path_stability <- path_stability[order(-path_stability$pi), ]

  return(path_stability = path_stability)
}

stability_scores <- stability(X = df_train, response = "Diagnosis")
