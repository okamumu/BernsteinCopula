#' @import mapfit
#' @title Phase-type marginal distribution constructor.
#'
#' @description Creates a MarginalDistribution object for a phase-type distribution.
#'
#' @param alpha A numeric vector of length k, representing the initial probabilities of the k states.
#' @param rate A numeric vector of length k, representing the rates of the k states.
#'
#' @return An object of class MarginalDistribution with appropriate pdf, cdf, and ccdf methods.
#'
#' @examples
#' fx <- phdist(alpha = c(0.5, 0.5), rate = c(1, 2))
#' fx$pdf(1)
#' fx$cdf(1)
#' fx$ccdf(1)
#'
#' @export
phdist <- function(alpha = c(0.5, 0.5), rate = c(1.0, 1.0)) {
  MarginalDistribution$new(
    name = "phasetype",
    pdf = function(x, params) dphase(x, ph = cf1(alpha = params$alpha, rate = params$rate)),
    cdf = function(x, params) pphase(x, ph = cf1(alpha = params$alpha, rate = params$rate), lower.tail = TRUE),
    ccdf = function(x, params) pphase(x, ph = cf1(alpha = params$alpha, rate = params$rate), lower.tail = FALSE),
    params = list(alpha = alpha, rate = rate),
    emstep_func = function(x, w1, w2, params) {
        data <- list()
        weights <- numeric()
        for (i in seq_along(x)) {
            xi <- x[i]
            data <- c(data, list(c(0, xi), xi, c(xi, Inf)))
            weights <- c(weights, w1[i], 1, w2[i])
        }
        mapfit::Remstep_cf1_interval(data, weights, params)
    }
  )
}

