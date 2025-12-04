build_paths <- function(X, response, K = 5, epsilon = .001, delta = 2, L = 25, model_type = NULL) {

  # --- Auto-detect model type ---
  detect_model_type <- function(data, response_name) {
    y <- data[[response_name]]

    if (is.factor(y) || length(unique(y)) == 2) return("logistic")
    if (is.numeric(y)) return("linear")
    stop("Cannot detect model type automatically.")
  }

  if (is.null(model_type)) {
    model_type <- detect_model_type(X, response)
    cat("Detected model type:", model_type, "\n")
  }

  predictors <- colnames(X)
  predictors <- setdiff(predictors, response)  # exclude response

  # --- Model fitting helper ---
  fit_model <- function(vars) {
    formula_str <- paste(response, "~", ifelse(length(vars) == 0, "1", paste(vars, collapse = "+")))
    if (model_type == "linear") {
      return(lm(as.formula(formula_str), data = X))
    } else {
      return(glm(as.formula(formula_str), data = X, family = binomial))
    }
  }

  # --- Initialize storage ---
  frontiers <- list()
  aic_by_model <- list()

  # --- Level 1: single-variable models ---
  start_fit <- vector("list", length(predictors))
  start_aic <- rep(NA, length(predictors))


  for (i in seq_along(predictors)) {
    suppressWarnings({
      start_fit[[i]] <- try(fit_model(predictors[i]), silent = TRUE)
    })
    if (!inherits(start_fit[[i]], "try-error")) {
      start_aic[i] <- AIC(start_fit[[i]])
    } else {
      start_fit[[i]] <- NA
      start_aic[i] <- NA
    }
  }

  # Save first frontier
  frontiers[[1]] <- lapply(seq_along(predictors), function(i) list(vars = predictors[i], fit = start_fit[[i]], aic = start_aic[i]))
  aic_by_model[[1]] <- start_aic

  # --- Loop over higher dimensions ---
  current_frontier <- frontiers[[1]]

  for (k in 2:K) {

    all_children <- list()
    all_children_aic <- numeric(0)

    for (parent in current_frontier) {
      parent_vars <- parent$vars
      parent_aic  <- parent$aic

      remaining <- setdiff(predictors, parent_vars)
      child_list <- list()
      child_aics <- numeric(0)

      for (v in remaining) {
        new_vars <- c(parent_vars, v)
        fit <- try(fit_model(new_vars), silent = TRUE)
        if (!inherits(fit, "try-error")) {
          aic <- AIC(fit)
          child_list <- append(child_list, list(list(vars = new_vars, fit = fit, aic = aic)))
          child_aics <- c(child_aics, aic)
        }
      }

      # Option A filtering
      if (length(child_aics) > 0) {
        best_child_aic <- min(child_aics, na.rm = TRUE)
        keep_idx <- which((child_aics - best_child_aic <= delta) & (parent_aic - child_aics > epsilon))
        if (length(keep_idx) > 0) {
          all_children <- append(all_children, child_list[keep_idx])
          all_children_aic <- c(all_children_aic, child_aics[keep_idx])
        }
      }
    }

    # Deduplicate by variable set
    if (length(all_children) > 0) {
      keys <- sapply(all_children, function(m) paste(sort(m$vars), collapse = "_"))
      all_children <- all_children[!duplicated(keys)]
      all_children_aic <- sapply(all_children, function(m) m$aic)

      # Keep top L by AIC if too many
      if (length(all_children) > L) {
        idx <- order(all_children_aic)[1:L]
        all_children <- all_children[idx]
      }
    }

    # Save this level
    frontiers[[k]] <- all_children
    aic_by_model[[k]] <- sapply(all_children, function(m) m$aic)
    current_frontier <- all_children
  }

  meta <- list(Dimensions = K, Epsilon = epsilon, Delta = delta, Limit_For_Retention = L, Model_Type = model_type)

  return(list(
    path_forest = list(
      frontiers = frontiers,
      aic_by_model = aic_by_model,
      meta = meta
    )
  ))
}

library(dplyr)

data("BreastCancer", package = "mlbench")
df <- BreastCancer

# Clean and preprocess
df <- df %>%
  select(-Id) %>%                         # Drop ID column
  filter(complete.cases(.)) %>%          # Remove rows with missing values
  mutate(across(-Class, as.numeric)) %>% # Convert predictors to numeric
  mutate(Class = factor(Class))          # Ensure response is factor

# Rename response to "Diagnosis" for clarity
names(df)[names(df) == "Class"] <- "Diagnosis"

test <- build_paths(X = df, response = "Diagnosis", K = 5, epsilon = .001, delta = 10, L = 25, model_type = NULL)

##############################################################################################

test_stability <- stability()

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

  return(list(path_stability = path_stability, resampling = resampling))
}






