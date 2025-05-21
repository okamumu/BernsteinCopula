#' Empirical Marginal Distribution Constructor
#'
#' Creates a \code{MarginalDistribution} object for an empirical distribution
#' based on the provided data.
#'
#' @param data Numeric vector. Data points for the empirical distribution.
#' @param ecdf Function. Empirical CDF function (default is \code{ecdf(data)}).
#' @param density Function. Density function (default is \code{density(data)}).
#'
#' @return A \code{MarginalDistribution} R6 object with fields \code{params}, \code{pdf}, \code{cdf}, \code{ccdf}, and \code{emstep_func}.
#'
#' @export
empdist <- function(data, ecdf = NULL, density = NULL) {
  if (is.null(ecdf)) {
    ecdf <- ecdf(data)
  }
  if (is.null(density)) {
    density <- density(data)
  }
  MarginalDistribution$new(
    name = "empirical",
    pdf = function(x, params) approx(params$density$x, params$density$y, xout = x, rule = 2)$y,
    cdf = function(x, params) params$ecdf(x),
    ccdf = function(x, params) 1 - params$ecdf(x),
    params = list(ecdf = ecdf, density = density),
    emstep_func = NULL,
    emstep_func_group = NULL
  )
}
