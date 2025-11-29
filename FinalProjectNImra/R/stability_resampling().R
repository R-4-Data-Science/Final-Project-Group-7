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
stability_resampling <- function(data, response, B = 50, K = 5,
                                 epsilon = 1e-6, delta = 2, L = 25,
                                 model_type = NULL) {

  # 🔍 Automatically detect model type if not provided
  if (is.null(model_type)) {
    y <- data[[response]]
    if (is.factor(y) || is.logical(y) || length(unique(y)) == 2) {
      model_type <- "logistic"
    } else {
      model_type <- "linear"
    }
    cat("Detected model type:", model_type, "\n")
  }

  predictors <- setdiff(names(data), response)
  stability_matrix <- matrix(0, nrow = B, ncol = length(predictors))
  colnames(stability_matrix) <- predictors

  for (b in 1:B) {
    cat("\nBootstrap sample", b, "of", B, "\n")
    boot_indices <- sample(1:nrow(data), replace = TRUE)
    boot_data <- data[boot_indices, ]

    models_b <- multi_path_forward(boot_data, response, K, epsilon, delta, L, model_type)
    vars_selected <- unique(unlist(lapply(models_b, function(m) m$variables)))
    stability_matrix[b, vars_selected] <- 1
  }

  stability_scores <- colMeans(stability_matrix)
  stability_df <- tibble::tibble(
    Variable = names(stability_scores),
    Stability = stability_scores
  )

  stability_df <- stability_df[order(-stability_df$Stability), ]
  return(stability_df)
}
