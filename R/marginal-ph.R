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
    pdf = function(x, params) {
      out <- rep(0, length(x))
      valid <- x > 0 & is.finite(x)
      out[valid] <- dphase(x[valid], ph = cf1(alpha = params$alpha, rate = params$rate))
      out
    },
    cdf = function(x, params) {
      out <- numeric(length(x))
      out[x > 0 & is.finite(x)] <-
        pphase(x[x > 0 & is.finite(x)],
              ph = cf1(alpha = params$alpha, rate = params$rate),
              lower.tail = TRUE)
      out[x == Inf] <- 1
      out
    },
    ccdf = function(x, params) {
      out <- numeric(length(x)) + 1
      out[x > 0 & is.finite(x)] <-
        pphase(x[x > 0 & is.finite(x)],
              ph = cf1(alpha = params$alpha, rate = params$rate),
              lower.tail = FALSE)
      out[x == Inf] <- 0
      out
    },
    params = list(alpha = alpha, rate = rate),
    emstep_func = function(x, w1, w2, params) {
      res <- emstep_ph_build_interval_data(x, w1, w2)
      b <- res$weights != 0.0
      mapfit::Remstep_cf1_interval(res$data[b], res$weights[b], params)
    },
    emstep_func_group = function(x, w1, w2, params) {
      res <- emstep_ph_build_group_data(x, w1, w2)
      b <- res$weights != 0.0
      mapfit::Remstep_cf1_interval(res$data[b], res$weights[b], params)
    }
  )
}

