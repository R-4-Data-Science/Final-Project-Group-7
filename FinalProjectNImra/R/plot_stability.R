#' Plot Variable Stability
#'
#' Visualizes stability scores as a bar plot.
#'
#' @param stability_df Tibble from stability_resampling().
#' @export
plot_stability <- function(stability_df) {
  ggplot2::ggplot(stability_df, ggplot2::aes(x = reorder(Variable, Stability), y = Stability)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Variable Stability Scores",
                  x = "Predictor",
                  y = "Stability (Proportion of Selection)") +
    ggplot2::theme_minimal(base_size = 14)
}

