#' Exponential Marginal Distribution Constructor
#'
#' Creates a \code{MarginalDistribution} object for an exponential distribution
#' with a specified rate parameter.
#'
#' @param rate Numeric. Rate parameter (lambda) of the exponential distribution (default is 1).
#'
#' @return An object of class \code{MarginalDistribution} with appropriate \code{pdf}, \code{cdf}, and \code{ccdf} methods.
#'
#' @examples
#' fx <- exponentialdist(rate = 2)
#' fx$pdf(1)
#' fx$cdf(1)
#' fx$ccdf(1)
#'
#' @export
exponentialdist <- function(rate = 1) {
  MarginalDistribution$new(
    name = "exponential",
    pdf = function(x, params) dexp(x, rate = params$rate),
    cdf = function(x, params) pexp(x, rate = params$rate, lower.tail = TRUE),
    ccdf = function(x, params) pexp(x, rate = params$rate, lower.tail = FALSE),
    params = list(rate = rate),
    emstep_func = emstep_exponential,
    emstep_func_group = emstep_exponential_group
  )
}
