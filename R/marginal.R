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
      name = NULL,
      params = list(),
      pdf_func = NULL,
      cdf_func = NULL,
      ccdf_func = NULL,

      #' @description
      #' Initialize a MarginalDistribution object.
      #' @param name Character. Name of the distribution (e.g., "normal").
      #' @param pdf Function. Custom PDF function.
      #' @param cdf Function. Custom CDF function.
      #' @param ccdf Function. Custom CCDF function.
      #' @param params Named list. Parameters for the distribution.
      #' @return A MarginalDistribution object.
      initialize = function(name = NULL, pdf = NULL, cdf = NULL, ccdf = NULL, params = list()) {
        self$name <- name
        self$params <- params
        self$pdf_func <- pdf
        self$cdf_func <- cdf
        self$ccdf_func <- ccdf
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
#' @return An object of class \code{MarginalDistribution} with appropriate \code{pdf}, \code{cdf}, and \code{ccdf} methods.
#'
#' @examples
#' fx <- normaldist(mean = 0, sd = 1)
#' fx$pdf(0)
#' fx$cdf(0)
#' fx$ccdf(0)
#'
#' @export
normaldist <- function(mean = 0, sd = 1) {
  MarginalDistribution$new(
    name = "normal",
    pdf = function(x, params) dnorm(x, mean = params$mean, sd = params$sd),
    cdf = function(x, params) pnorm(x, mean = params$mean, sd = params$sd, lower.tail = TRUE),
    ccdf = function(x, params) pnorm(x, mean = params$mean, sd = params$sd, lower.tail = FALSE),
    params = list(mean = mean, sd = sd)
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
    params = list(rate = rate)
  )
}
