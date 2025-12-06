#' @title Stability Resampling
#'
#' @description This function estimates the stability of variable selection using repeated bootstrap method.
#'
#' @param X Data frame of predictors and response.
#' @param response is the name of the response variable
#' @param B Number of bootstrap samples.
#' @param K Maximum number of forward steps (default 5).
#' @param epsilon Minimum AIC improvement to continue from a parent (default 1e-6)
#' @param delta AIC tie threshold; keep children within `best + delta` (default 2)
#' @param L Maximum number of models to keep per step (default 25)
#' @param model_type Model type is specified ("linear" or "logistic"). Automatically detected if not provided ("NULL").
#' @return A dataframe 'path_stability' with two columns
#' \describe{
#'  \item{variable}{predictor names}
#'  \item{pi}{The stabilty score for each variable}
#' }
#' @export
#' @examples
#' # Stability function for a linear response
#' Linear_scores <- stability(X = mtcars, response = "mpg", K = 3, delta = 2)
#' # Stability function for a logistic response
#' mtcars$am <- factor(mtcars$am) #convert to factor
#' Logistic_scores <- stability(X = mtcars, response = "am", K = 3, delta = 2)
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
