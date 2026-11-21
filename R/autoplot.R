#' Plot a balance check
#'
#' Plot the output of [balance_check()]. `type = "smd"` shows the
#' standardized mean difference of each variable in each split, with dashed
#' lines at -0.1 and 0.1. For a three-way split, the validation and testing
#' sets have different shapes. `type = "energy"` shows the energy distance of
#' each split.
#'
#' @param object A `balance_check` object.
#' @param type `"smd"` or `"energy"`.
#' @param ... Not used.
#'
#' @return A ggplot object.
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' library(ggplot2)
#' set.seed(1)
#' folds <- balanced_vfold_cv(mtcars, v = 4, balance = c(mpg, hp, wt))
#' autoplot(balance_check(folds))
#' autoplot(balance_check(folds), type = "energy")
#' @exportS3Method ggplot2::autoplot
autoplot.balance_check <- function(object, type = c("smd", "energy"), ...) {
  rlang::check_installed("ggplot2", reason = "to plot a balance check.")
  type <- rlang::arg_match(type)
  df <- tibble::as_tibble(unclass_balance(object))
  has_id2 <- "id2" %in% names(df)

  if (type == "smd") {
    df <- df[df$statistic == "smd", , drop = FALSE]
    df$label <- ifelse(
      is.na(df$level),
      df$variable,
      paste0(df$variable, ": ", df$level)
    )
    df$label <- factor(df$label, levels = rev(unique(df$label)))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$label)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey50") +
      ggplot2::geom_vline(
        xintercept = c(-0.1, 0.1),
        linetype = "dashed",
        colour = "grey50"
      ) +
      ggplot2::labs(
        x = "Standardized mean difference (assessment - analysis)",
        y = NULL
      )
    # A three-way split has two comparisons, told apart by shape.
    if (!has_id2 && setequal(df$id, c("Validation", "Testing"))) {
      df$id <- factor(df$id, levels = c("Validation", "Testing"))
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(shape = .data$id),
          data = df,
          alpha = 0.6
        ) +
        ggplot2::labs(shape = NULL)
    } else {
      p <- p + ggplot2::geom_point(alpha = 0.6)
    }
    return(p)
  }

  df <- df[df$statistic == "energy", , drop = FALSE]
  if (has_id2) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$id2, y = .data$value, colour = .data$id)
    ) +
      ggplot2::geom_point() +
      ggplot2::labs(x = NULL, y = "Energy distance", colour = NULL)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$id, y = .data$value)) +
      ggplot2::geom_point() +
      ggplot2::labs(x = NULL, y = "Energy distance")
  }
  p + ggplot2::expand_limits(y = 0)
}

unclass_balance <- function(x) {
  class(x) <- setdiff(class(x), "balance_check")
  attr(x, "id_cols") <- NULL
  x
}
