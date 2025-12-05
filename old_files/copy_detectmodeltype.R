#' Detect Model Type
#'
#' Automatically detects whether the response variable should be modeled as
#' a linear or logistic regression based on its type.
#'
#' @param data A data frame.
#' @param response The name of the response variable.
#' @return A character string: "linear" or "logistic".
#' @export
detect_model_type <- function(data, response) {
  y <- data[[response]]
  
  # Handle factors
  if (is.factor(y)) {
    y_levels <- length(levels(y))
    if (y_levels == 2) return("logistic")
    else return("linear")
  }
  
  # Handle numeric
  if (is.numeric(y)) {
    unique_vals <- length(unique(y))
    if (unique_vals == 2) return("logistic")
    else return("linear")
  }
  
  # Handle logical
  if (is.logical(y)) return("logistic")
  
  # Default fallback
  return("linear")
}
