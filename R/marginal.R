#' @title MarginalDistribution R6 Class
#'
#' @description
#' R6 class for specifying univariate marginal distributions, including standard
#' ones (e.g., normal) and custom distributions.
#'
#' @docType class
#' @format An R6 class object with the following fields and methods:
#' \describe{
#'   \item{name}{Character. Name of the distribution.}
#'   \item{params}{Named list. Parameters passed to pdf/cdf.}
#'   \item{pdf, cdf, ccdf}{Function objects of the distribution.}
#'   \item{eval_pdf(x)}{Evaluate the PDF at point(s) x.}
#'   \item{eval_cdf(x)}{Evaluate the CDF at point(s) x.}
#'   \item{eval_ccdf(x)}{Evaluate the complementary CDF at point(s) x.}
#'   \item{print()}{Display the distribution name and parameters.}
#' }
#'
#' @export
MarginalDistribution <-
  R6::R6Class("MarginalDistribution",
    public = list(
      #' @field name Character. Name of the distribution.
      #' @field params Named list. Parameters passed to pdf/cdf.
      #' @field pdf_func Function. PDF function.
      #' @field cdf_func Function. CDF function.
      #' @field ccdf_func Function. CCDF function.
      #' @field emstep_func Function. EM step function for parameter estimation.
      name = NULL,
      params = list(),
      pdf_func = NULL,
      cdf_func = NULL,
      ccdf_func = NULL,
      emstep_func = NULL,

      #' @description
      #' Initialize a MarginalDistribution object.
      #' @param name Character. Name of the distribution (e.g., "normal").
      #' @param pdf Function. Custom PDF function.
      #' @param cdf Function. Custom CDF function.
      #' @param ccdf Function. Custom CCDF function.
      #' @param params Named list. Parameters for the distribution.
      #' @param emstep_func Function. EM step function for parameter estimation.
      #' @return A MarginalDistribution object.
      initialize = function(name = NULL, pdf = NULL, cdf = NULL, ccdf = NULL, params = list(), 
                            emstep_func = NULL) {
        self$name <- name
        self$params <- params
        self$pdf_func <- pdf
        self$cdf_func <- cdf
        self$ccdf_func <- ccdf
        self$emstep_func <- emstep_func
      },

      #' @description
      #' Evaluate the PDF at point(s) x.
      #' @param x Numeric vector. Points at which to evaluate the PDF.
      #' @return Numeric vector. PDF
      pdf = function(x) self$pdf_func(x, self$params),

      #' @description
      #' Evaluate the CDF at point(s) x.
      #' @param x Numeric vector. Points at which to evaluate the CDF.
      #' @return Numeric vector. CDF
      cdf = function(x) self$cdf_func(x, self$params),

      #' @description
      #' Evaluate the CCDF at point(s) x.
      #' @param x Numeric vector. Points at which to evaluate the CCDF.
      #' @return Numeric vector. CCDF
      ccdf = function(x) self$ccdf_func(x, self$params),

      #' @description
      #' Exec the EM step function.
      #' @param x Numeric vector. Data points.
      #' @param w1 Numeric vector. Weights for the data points.
      #' @param w2 Numeric vector. Weights for the data points.
      #' @return Named list. Updated parameters.
      emstep = function(x, w1, w2) {
        if (is.null(self$emstep_func)) {
          stop("EM step function is not defined for this distribution.")
        }
        self$emstep_func(x, w1, w2, self$params)
      },

      #' @description
      #' Print the distribution name and parameters.
      #' @param ... Additional arguments (not used).
      #' @return NULL
      print = function(...) {
        cat(sprintf("<MarginalDistribution: %s>\n", self$name))
        print(self$params)
      }
    )
  )

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
    emstep_func = emstep_normal
  )
}

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
    emstep_func = NULL
  )
}

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
    emstep_func = emstep_exponential
  )
}
