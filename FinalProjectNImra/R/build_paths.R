#' Multi-Path Forward Selection (auto-detects model type)
#'
#' Implements forward selection with multiple near-optimal AIC paths for
#' linear (Gaussian) or logistic (binomial) regression. If `model_type = NULL`,
#' the type is auto-detected from the response:
#' - Binary factor / logical / numeric in {0,1} → `"logistic"`
#' - Otherwise → `"linear"`
#'
#' @param data Data frame.
#' @param response Response variable name (character).
#' @param K Maximum number of forward steps (default 5).
#' @param epsilon Minimum AIC improvement to continue from a parent (default 1e-6).
#' @param delta AIC tie threshold; keep children within `best + delta` (default 2).
#' @param L Maximum number of models to keep per step (default 25).
#' @param model_type "linear", "logistic", or `NULL` for auto-detect (default `NULL`).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{path_forest}: list of frontiers (models kept at each step)
#'   \item \code{aic_by_model}: named list mapping "var1,var2,..." → AIC
#'   \item \code{meta}: list of meta settings (K, epsilon, delta, L, model_type)
#' }
#' @export
#'
#' @examples
#' # Linear (auto-detected): response is numeric with many values
#' mp_lin <- multi_path_forward(mtcars, response = "mpg", K = 3, delta = 2)
#' str(mp_lin$path_forest, max.level = 1)
#'
#' # Logistic (auto-detected): response has 2 levels
#' mt <- within(mtcars, { am <- factor(am) })
#' mp_log <- multi_path_forward(mtcars, response = "am", K = 3, delta = 2)
#' str(mp_log$path_forest, max.level = 1)
#'
multi_path_forward <- function(data, response,
                               K = 5, epsilon = 1e-6, delta = 2, L = 25,
                               model_type = NULL) {

  # ---- Local helper: detect model type -------------------------------------
  detect_model_type <- function(y) {
    if (is.factor(y) && nlevels(y) == 2) return("logistic")
    if (is.logical(y)) return("logistic")
    if (is.numeric(y)) {
      uy <- unique(stats::na.omit(y))
      if (length(uy) <= 2 && all(uy %in% c(0, 1))) return("logistic")
    }
    "linear"
  }

  # ---- Auto-detect model type ----------------------------------------------
  if (is.null(model_type)) {
    model_type <- detect_model_type(data[[response]])
    cat("Detected model type:", model_type, "\n")
  } else {
    model_type <- match.arg(tolower(model_type), c("linear", "logistic"))
  }

  predictors <- setdiff(names(data), response)

  # ---- Model fitting helper ------------------------------------------------
  fit_model <- function(vars) {
    formula_str <- paste(
      response, "~",
      ifelse(length(vars) == 0, "1", paste(vars, collapse = "+"))
    )
    f <- stats::as.formula(formula_str)

    if (model_type == "linear") {
      stats::lm(f, data = data)
    } else {
      stats::glm(f, data = data, family = stats::binomial())
    }
  }

  # ---- Initialize with intercept-only model --------------------------------
  start_fit <- fit_model(character())
  start_aic <- stats::AIC(start_fit)

  models <- list(
    list(
      variables = character(),
      AIC = start_aic,
      fit = start_fit
    )
  )

  frontiers <- list()
  aic_by_model <- list()

  step <- 1L

  repeat {
    new_models <- list()
    model_hash <- character()

    for (m in models) {
      current_vars <- m$variables
      remaining_vars <- setdiff(predictors, current_vars)
      if (length(remaining_vars) == 0L) next

      candidate_models <- list()

      for (v in remaining_vars) {
        new_vars <- c(current_vars, v)
        key <- paste(sort(new_vars), collapse = ",")

        if (key %in% model_hash) next

        fit <- try(
          suppressWarnings(fit_model(new_vars)),
          silent = TRUE
        )
        if (inherits(fit, "try-error")) next

        aic_val <- stats::AIC(fit)
        if (!is.finite(aic_val)) next

        candidate_models[[key]] <- list(
          variables = new_vars,
          AIC = aic_val,
          fit = fit
        )
        model_hash <- c(model_hash, key)
      }

      if (length(candidate_models) == 0L) next

      aic_vals <- vapply(candidate_models, `[[`, numeric(1), "AIC")
      best_aic <- min(aic_vals)

      # ---- Keep near-optimal children ---------------------------------------
      selected <- candidate_models[aic_vals <= best_aic + delta]

      # ---- Progress only if parent improves ---------------------------------
      if (best_aic < m$AIC - epsilon) {
        new_models <- c(new_models, selected)
      }
    }

    if (length(new_models) == 0L || step >= K) break

    # ---- Sort & truncate ----------------------------------------------------
    aics <- vapply(new_models, `[[`, numeric(1), "AIC")
    ord <- order(aics)
    new_models <- new_models[ord]
    new_models <- new_models[seq_len(min(L, length(new_models)))]

    frontiers[[step]] <- new_models

    # ---- Save AIC dictionary ------------------------------------------------
    for (m in new_models) {
      key <- paste(sort(m$variables), collapse = ",")
      aic_by_model[[key]] <- m$AIC
    }

    models <- new_models
    step <- step + 1L
  }

  list(
    path_forest = frontiers,
    aic_by_model = aic_by_model,
    meta = list(
      K = K,
      epsilon = epsilon,
      delta = delta,
      L = L,
      model_type = model_type
    )
  )
}

