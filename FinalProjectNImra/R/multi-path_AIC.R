#' Full Multi-Path AIC Procedure
#'
#' Combines multi-path model search, stability estimation, and plausible model selection
#' into a complete end-to-end workflow.
#'
#' @param data A data frame containing predictors and response.
#' @param response Character name of the response variable.
#' @param K,epsilon,delta,L Parameters controlling the path search.
#' @param B Number of bootstrap samples for stability estimation.
#' @param delta_AIC,tau Thresholds for plausible model selection.
#' @param model_type Optional model type ("linear" or "logistic").
#' If not provided, detected automatically.
#'
#' @return A list containing:
#'   \item{stability}{Variable stability tibble.}
#'   \item{plausible_models}{List of plausible models.}
#' @export

multi_path_AIC_procedure <- function(data, response,
                                     K = 5, epsilon = 1e-6, delta = 2, L = 25,
                                     B = 50, delta_AIC = 2, tau = 0.6,
                                     model_type = NULL) {

  # Automatically detect model type if not specified
  if (is.null(model_type)) {
    model_type <- detect_model_type(data, response)
    cat("Detected model type:", model_type, "\n")
  }

  cat("=== STEP 1: Multi-path search on full data ===\n")
  full_data_models <- multi_path_forward(data, response, K, epsilon, delta, L, model_type)

  cat("\n=== STEP 2: Stability estimation with resampling ===\n")
  stability_results <- stability_resampling(data, response, B, K, epsilon, delta, L, model_type)
  print(stability_results)

  cat("\n=== STEP 3: Selecting plausible, stable models ===\n")
  plausible_models <- select_plausible_models(
    full_data_models = full_data_models,
    stability_df = stability_results,
    delta_AIC = delta_AIC,
    tau = tau
  )

  cat("\n=== DONE ===\n")

  return(list(
    stability = stability_results,
    plausible_models = plausible_models
  ))
}
