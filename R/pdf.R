#' Bernstein Copula PDF
#'
#' Computes the probability density function of the Bernstein copula at a given point (u, v).
#'
#' @param u A numeric value in [0, 1], the first copula argument.
#' @param v A numeric value in [0, 1], the second copula argument.
#' @param R A numeric matrix of size m x n representing the copula weights.
#'
#' @return A numeric value giving the copula density at (u, v).
#' @export
#'
#' @examples
#' R <- matrix(1/4, 2, 2)
#' bernstein_copula_pdf(0.5, 0.5, R)
#'
#' @useDynLib BernsteinCopula, .registration = TRUE
#' @importFrom Rcpp evalCpp
bernstein_copula_pdf <- function(u, v, R) {
  .Call(`_BernsteinCopula_bernstein_copula_pdf`, u, v, R)
}
