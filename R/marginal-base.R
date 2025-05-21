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
      #' @field emstep_func_group Function. EM step function for parameter estimation in groups.
      name = NULL,
      params = list(),
      pdf_func = NULL,
      cdf_func = NULL,
      ccdf_func = NULL,
      emstep_func = NULL,
      emstep_func_group = NULL,

      #' @description
      #' Initialize a MarginalDistribution object.
      #' @param name Character. Name of the distribution (e.g., "normal").
      #' @param pdf Function. Custom PDF function.
      #' @param cdf Function. Custom CDF function.
      #' @param ccdf Function. Custom CCDF function.
      #' @param params Named list. Parameters for the distribution.
      #' @param emstep_func Function. EM step function for parameter estimation.
      #' @param emstep_func_group Function. EM step function for parameter estimation in groups.
      #' @return A MarginalDistribution object.
      initialize = function(name = NULL, pdf = NULL, cdf = NULL, ccdf = NULL, params = list(), 
                            emstep_func = NULL, emstep_func_group = NULL) {
        self$name <- name
        self$params <- params
        self$pdf_func <- pdf
        self$cdf_func <- cdf
        self$ccdf_func <- ccdf
        self$emstep_func <- emstep_func
        self$emstep_func_group <- emstep_func_group
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
