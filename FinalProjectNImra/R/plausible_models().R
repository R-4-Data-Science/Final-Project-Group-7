#' @title Plausible Model Selection with AIC and Stability Filtering
#'
#' @description
#' Selects plausible and stable models using AIC and average stability thresholds
#'
#' @param full_data_models List of models from multi_path_forward().
#' @param stability_df Tibble of stability scores from stability_resampling().
#' @param delta Numeric; maximum AIC difference from best model.(default 2)
#' @param tau Numeric; minimum average stability threshold.(default 0.6)
#' @return A List ("plausible_models") of plausible and stable models that contain:
#' \describe{
#'  \item{fit}{A fitted model}
#'  \item{aic}{AIC value for the model}
#'  \item{vars}{Predictor variable names within models}
#'  \item{avg_stability}{Average stability score of the predicitor models}
#' }
#' @export
plausible_models <- function(full_data_models, stability_df,
                             delta = 2, tau = 0.6) {

  # unlist models, predictors and aic values from full set of models generated in previous step,
  # combine into one master list
  unlist_models <- list()
  collect_models <- list()
  collect_aic <- list()
  collect_predictors <- list()

  unlist_models <- unlist(full_data_models, recursive = FALSE)

  for(b in 1:length(unlist_models)){

    collect_models[[b]] <- unlist_models[[b]]$fit
    collect_aic[[b]] <- unlist_models[[b]]$aic
    collect_predictors[[b]] <- unlist_models[[b]]$vars

  }

  # find minimum model aic
  min_aic <- min(unlist(collect_aic))

  # make list to save plausible models
  plausible_models <- list()

  for(b in 1:length(unlist_models)){
    if(unlist_models[[b]]$aic <= (min_aic + delta)){
      plausible_models <- append(plausible_models, list(unlist_models[[b]]))
    }
  }

  # Compute mean stability for each model
  stability_vec <- setNames(stability_scores$pi, stability_scores$variable)
  for (i in seq_along(plausible_models)) {
    vars <- plausible_models[[i]]$vars
    avg_stab <- mean(stability_vec[vars], na.rm = TRUE)
    plausible_models[[i]]$avg_stability <- avg_stab
  }

  # Filter by average stability threshold
  plausible_models <- Filter(function(m) m$avg_stability >= tau, plausible_models)

  return(plausible = plausible_models)
}
