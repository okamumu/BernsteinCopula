#' Bin counts for 2D data
#'
#' @description This function takes two numeric vectors and bins them into a 2D histogram.
#' It returns a matrix of counts for each bin.
#'
#' @param x Numeric vector. First variable to be binned.
#' @param y Numeric vector. Second variable to be binned.
#' @param xbreaks Numeric or character. Breaks for the x variable. Default is "Sturges".
#' @param ybreaks Numeric or character. Breaks for the y variable. Default is "Sturges".
#'
#' @return A list containing:
#' \item{xbreaks}{Breaks for the x variable.}
#' \item{ybreaks}{Breaks for the y variable.}
#' \item{counts}{A matrix of counts for each bin.}
#' The dimensions of the matrix correspond to the number of bins in x and y.
#' The row names correspond to the x bins and the column names correspond to the y bins.
#'
#' @examples
#' x <- rnorm(1000)
#' y <- rnorm(1000)
#' counts <- bincount2d(x, y)
#' print(counts)
#' @export
#' @importFrom graphics hist
bincount2d <- function(x, y, xbreaks = "Sturges", ybreaks = "Sturges") {
  stopifnot(length(x) == length(y))
  xbreaks <- hist(x, breaks = xbreaks, plot = FALSE)$breaks
  ybreaks <- hist(y, breaks = ybreaks, plot = FALSE)$breaks
  xbin <- cut(x, breaks = xbreaks, include.lowest = TRUE, right = FALSE)
  ybin <- cut(y, breaks = ybreaks, include.lowest = TRUE, right = FALSE)
  counts <- table(xbin, ybin)
  list(
    xbreaks = xbreaks,
    ybreaks = ybreaks,
    counts = as.matrix(counts)
  )
}