#' Joint PDF from Bernstein Copula and MarginalDistributions
#'
#' Computes the joint density \eqn{h(x, y) = c(F(x), G(y)) \cdot f(x) \cdot g(y)} using a Bernstein copula and marginal distributions.
#'
#' @param x Numeric scalar. Evaluation point for variable X.
#' @param y Numeric scalar. Evaluation point for variable Y.
#' @param R Matrix of size m x n. Copula weight matrix.
#' @param Fx A \code{MarginalDistribution} object for variable X.
#' @param Gy A \code{MarginalDistribution} object for variable Y.
#'
#' @return Numeric. Value of the joint density \code{h(x, y)}.
#'
#' @export
#'
#' @examples
#' Fx <- normaldist(mean = 0, sd = 1)
#' Gy <- normaldist(mean = 0, sd = 1)
#' R <- matrix(1/4, 2, 2)
#' bernstein_joint_pdf(0.5, 0.5, R, Fx, Gy)
bernstein_joint_pdf <- function(x, y, R, Fx, Gy) {
  u <- Fx$cdf(x)
  v <- Gy$cdf(y)
  copula_density <- bernstein_copula_pdf(u, v, R)
  Fx$pdf(x) * Gy$pdf(y) * copula_density
}

#' Joint PDF Grid from Bernstein Copula and MarginalDistributions
#'
#' Evaluates the joint density \eqn{h(x, y)} over a grid of values using a Bernstein copula and two \code{MarginalDistribution} objects.
#'
#' @param x Numeric vector. Grid of x values.
#' @param y Numeric vector. Grid of y values.
#' @param R Matrix of size m x n. Bernstein copula weight matrix.
#' @param Fx A \code{MarginalDistribution} object for variable X.
#' @param Gy A \code{MarginalDistribution} object for variable Y.
#'
#' @return A matrix of joint density values \code{h(x, y)} with dimensions \code{length(x) × length(y)}.
#'
#' @export
#'
#' @examples
#' x <- seq(-2, 2, length.out = 50)
#' y <- seq(-2, 2, length.out = 50)
#' Fx <- normaldist(mean = 0, sd = 1)
#' Gy <- normaldist(mean = 0, sd = 1)
#' R <- matrix(1/4, 2, 2)
#' Z <- bernstein_joint_pdf_grid(x, y, R, Fx, Gy)
#' plot_contour(x, y, Z)
bernstein_joint_pdf_grid <- function(x, y, R, Fx, Gy) {
  u_vals <- Fx$cdf(x)
  v_vals <- Gy$cdf(y)
  fx_vals <- Fx$pdf(x)
  gy_vals <- Gy$pdf(y)

  outer_result <- outer(seq_along(x), seq_along(y), Vectorize(function(i, j) {
    u <- u_vals[i]
    v <- v_vals[j]
    copula <- bernstein_copula_pdf(u, v, R)
    fx_vals[i] * gy_vals[j] * copula
  }))

  dimnames(outer_result) <- list(x = x, y = y)
  return(outer_result)
}
