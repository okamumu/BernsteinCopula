#' M-step for Bernstein Copula EM Algorithm (Dou et al.)
#'
#' Performs the M-step of the EM algorithm described in Dou et al. (2016),
#' estimating the Bernstein copula weight matrix by solving a constrained
#' optimization problem via iterative bisection.
#'
#' @param tau_bar A numeric matrix of pseudo-counts (typically from E-step).
#' @param maxiter Integer. Maximum number of iterations. Default is 1000.
#' @param tol Numeric. Convergence tolerance. Default is 1e-10.
#'
#' @return A normalized matrix \code{R} such that \code{sum(R) = 1} and
#'         \code{R[i,j] ≈ tau_bar[i,j] / (mu[i] + lambda[j])}.
#' @export
#'
#' @examples
#' set.seed(1)
#' tau <- matrix(runif(16), 4, 4)
#' R <- dou_em_mstep(tau)
#' sum(R)  # should be approximately 1
dou_em_mstep <- function(tau_bar, maxiter = 1000, tol = 1e-10) {
  .Call(`_BernsteinCopula_dou_mstep`, tau_bar, maxiter, tol)
}
