#' Contour Plot from Grid Data Using ggplot2
#'
#' Creates a contour or filled contour plot using ggplot2 from grid data (like base R's contour).
#'
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param z Matrix of z values. Dimensions must match length(x) × length(y).
#' @param filled Logical. If TRUE, use filled contour (`geom_contour_filled`). Default FALSE.
#' @param title Character. Plot title. Default: "Contour Plot".
#' @param xlabel Character. Label for x axis. Default: "x".
#' @param ylabel Character. Label for y axis. Default: "y".
#' @param xlim Numeric vector of length 2, giving x axis limits. Default NULL (auto).
#' @param ylim Numeric vector of length 2, giving y axis limits. Default NULL (auto).
#' @param ... Additional arguments passed to the ggplot2 layer.
#'
#' @return A ggplot object with the contour plot.
#' @export
#'
#' @examples
#' x <- seq(-3, 3, length.out = 100)
#' y <- seq(-3, 3, length.out = 100)
#' z <- outer(x, y, function(x, y) dnorm(x) * dnorm(y))
#' plot_contour(x, y, z, xlim = c(-1, 1), ylim = c(-2, 2),
#'                 title = "Zoomed PDF", xlabel = "X", ylabel = "Y")
plot_contour <- function(x, y, z,
                            filled = FALSE,
                            title = "Contour Plot",
                            xlabel = "x",
                            ylabel = "y",
                            xlim = NULL,
                            ylim = NULL,
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

  p <- p + ggplot2::labs(title = title, x = xlabel, y = ylabel) +
    ggplot2::theme_minimal()

  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }

  return(p)
}
