#' Sinkhorn Scaling for Bernstein Copula Matrix
#'
#' Scales a non-negative matrix \code{W} so that its row and column sums match the Bernstein copula normalization constraints.
#' After convergence, the scaled matrix will have row sums equal to \code{1 / nrow(W)} and column sums equal to \code{1 / ncol(W)}.
#'
#' @param W A numeric matrix with non-negative entries.
#' @param maxiter Integer. Maximum number of iterations. Default is 1000.
#' @param tol Numeric. Convergence tolerance. Default is 1e-10.
#'
#' @return A scaled matrix of the same dimension as \code{W}, normalized to meet the marginal constraints.
#' @export
#'
#' @examples
#' W <- matrix(runif(16), 4, 4)
#' R <- sinkhorn_scaling(W)
#' round(rowSums(R), 6)  # Should be close to rep(1/4, 4)
#' round(colSums(R), 6)  # Should be close to rep(1/4, 4)
sinkhorn_scaling <- function(W, maxiter = 1000, tol = 1e-10) {
  .Call(`_BernsteinCopula_sinkhorn_scaling`, W, maxiter, tol)
}
