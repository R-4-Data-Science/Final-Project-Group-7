#' Multi-Path Forward Selection
#'
#' Implements forward selection with multiple near-optimal AIC paths.
#'
#' @param data Data frame.
#' @param response Response variable name.
#' @param K Maximum steps.
#' @param epsilon Minimum improvement threshold.
#' @param delta AIC tie threshold.
#' @param L Maximum number of models to keep.
#' @param model_type "linear" or "logistic".
#' @return A list of selected models and their AIC values.
#' @export
multi_path_forward <- function(data, response,
                               K = 5, epsilon = 1e-6, delta = 2, L = 25,
                               model_type = "linear") {
  predictors <- setdiff(names(data), response)

  fit_model <- function(vars) {
    formula_str <- paste(response, "~", ifelse(length(vars) == 0, "1", paste(vars, collapse = "+")))
    if (model_type == "linear") {
      return(lm(as.formula(formula_str), data = data))
    } else if (model_type == "logistic") {
      return(glm(as.formula(formula_str), data = data, family = binomial))
    } else {
      stop("Invalid model_type. Use 'linear' or 'logistic'.")
    }
  }

  start_fit <- fit_model(c())
  models <- list(list(variables = c(), AIC = AIC(start_fit)))

  step <- 1
  repeat {
    cat("\nStep", step, "\n")
    new_models <- list()

    for (m in models) {
      current_vars <- m$variables
      remaining_vars <- setdiff(predictors, current_vars)
      if (length(remaining_vars) == 0) next

      candidate_models <- list()
      for (v in remaining_vars) {
        new_vars <- c(current_vars, v)
        suppressWarnings({
          fit <- fit_model(new_vars)
        })
        aic_val <- AIC(fit)
        candidate_models[[v]] <- list(variables = new_vars, AIC = aic_val)
      }

      aic_values <- sapply(candidate_models, function(x) x$AIC)
      best_aic <- min(aic_values)
      near_tie_models <- candidate_models[aic_values <= best_aic + delta]

      if (best_aic < m$AIC - epsilon) {
        new_models <- c(new_models, near_tie_models)
      }
    }

    if (length(new_models) == 0 || step >= K) {
      cat("No further AIC improvement. Stopping at step", step, "\n")
      break
    }

    models <- new_models[1:min(L, length(new_models))]
    cat("Number of models kept:", length(models), "\n")
    for (i in seq_along(models)) {
      cat("Model", i, ":", paste(models[[i]]$variables, collapse = ", "),
          "| AIC =", round(models[[i]]$AIC, 3), "\n")
    }
    step <- step + 1
  }

  return(models)
}
