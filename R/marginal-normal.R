#' Normal Marginal Distribution Constructor
#'
#' Creates a \code{MarginalDistribution} object for a normal (Gaussian) distribution,
#' with specified mean and standard deviation.
#'
#' @param mean Numeric. Mean of the normal distribution (default is 0).
#' @param sd Numeric. Standard deviation of the normal distribution (default is 1).
#'
#' @return A \code{MarginalDistribution} R6 object with fields \code{params}, \code{pdf}, \code{cdf}, \code{ccdf}, and \code{emstep_func}.
#'
#' @export
normaldist <- function(mean = 0, sd = 1) {
  MarginalDistribution$new(
    name = "normal",
    pdf = function(x, params) dnorm(x, mean = params$mean, sd = params$sd),
    cdf = function(x, params) pnorm(x, mean = params$mean, sd = params$sd, lower.tail = TRUE),
    ccdf = function(x, params) pnorm(x, mean = params$mean, sd = params$sd, lower.tail = FALSE),
    params = list(mean = mean, sd = sd),
    emstep_func = emstep_normal,
    emstep_func_group = emstep_normal_group
  )
}
