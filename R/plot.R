#' Contour Plot from Grid Data Using ggplot2
#'
#' Creates a contour or filled contour plot using ggplot2 from grid data (like base R's contour),
#' with optional overlay of points.
#'
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param z Matrix of z values. Dimensions must match length(x) × length(y).
#' @param filled Logical. If TRUE, use filled contour (`geom_contour_filled`). Default TRUE.
#' @param title Character. Plot title. Default: "Contour Plot".
#' @param xlabel Character. Label for x axis. Default: "x".
#' @param ylabel Character. Label for y axis. Default: "y".
#' @param xlim, ylim Axis limits. Default NULL.
#' @param points Optional data frame with columns `x` and `y` to overlay on the plot.
#' @param ... Additional arguments passed to the ggplot2 layer.
#'
#' @return A ggplot object.
#' @export
plot_contour <- function(x, y, z,
                         filled = TRUE,
                         title = "Contour Plot",
                         xlabel = "x",
                         ylabel = "y",
                         xlim = NULL,
                         ylim = NULL,
                         points = NULL,
                         ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required but not installed.")
  }

  if (!is.matrix(z) || length(x) != nrow(z) || length(y) != ncol(z)) {
    stop("z must be a matrix with dimensions length(x) × length(y)")
  }

  df <- expand.grid(x = x, y = y)
  df$z <- as.vector(z)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, z = z))
  if (filled) {
    p <- p + ggplot2::geom_contour_filled(...)
  } else {
    p <- p + ggplot2::geom_contour(...)
  }

  if (!is.null(points)) {
    if (!all(c("x", "y") %in% names(points))) {
      stop("points must be a data frame with columns 'x' and 'y'")
    }
    p <- p + ggplot2::geom_point(data = points,
                                 ggplot2::aes(x = x, y = y),
                                 inherit.aes = FALSE,
                                 alpha = 0.4, color = "black")
  }

  p <- p + ggplot2::labs(title = title, x = xlabel, y = ylabel) +
    ggplot2::theme_minimal()

  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }

  return(p)
}
