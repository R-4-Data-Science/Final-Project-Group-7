#' @title Confusion Metrics
#' @description Compute confusion matrix and standard classification metrics
#' for logistic regression models at cutoff = 0.5.
#' @param model A fitted logistic regression model (`glm` with family = binomial)
#' @param data Data frame used for prediction
#' @param threshold Numeric, classification cutoff (default = 0.5)
#' @return A \code{list} containing confusion matrix and key metrics:
#' accuracy, sensitivity, specificity, precision, F1 score, and DOR.
#' @export
#' @examples
#' fit <- glm(am ~ wt + hp, data = mtcars, family = binomial)
#' confusion_metrics(fit, mtcars)
confusion_metrics <- function(model, data, threshold = 0.5) {
  if (!inherits(model, "glm") || family(model)$family != "binomial") {
    stop("Model must be a logistic regression (binomial family).")
  }

  probs <- predict(model, newdata = data, type = "response")
  pred <- ifelse(probs > threshold, 1, 0)
  actual <- data[[all.vars(formula(model))[1]]]
  ref <- levels(actual)[1]         # reference level (0)
  actual <- as.integer(actual != ref)

  TP <- sum(pred == 1 & actual == 1)
  TN <- sum(pred == 0 & actual == 0)
  FP <- sum(pred == 1 & actual == 0)
  FN <- sum(pred == 0 & actual == 1)

  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  precision <- TP / (TP + FP)
  F1 <- (2 * precision * sensitivity) / (precision + sensitivity)
  DOR <- (TP / FP) / (FN / TN)

  confusion_matrix <- matrix(c(TN, FP, FN, TP), nrow = 2,
                             dimnames = list("Actual" = c("0", "1"),
                                             "Predicted" = c("0", "1")))

  return(list(
    confusion_matrix = confusion_matrix,
    metrics = list(
      Accuracy = accuracy,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Precision = precision,
      F1 = F1,
      DOR = DOR
    )
  ))
}
