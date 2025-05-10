#' Run EM Algorithm for Bernstein Copula Estimation (R wrapper)
#'
#' Wrapper for the C++ EM implementation. Estimates the copula weight matrix \code{R}
#' given pseudo-observations \code{u}, \code{v}, and an initial matrix \code{R_init}.
#'
#' @param u Numeric vector in [0, 1]. Pseudo-observations for variable X.
#' @param v Numeric vector in [0, 1]. Pseudo-observations for variable Y.
#' @param R_init Numeric matrix (m × n). Initial Bernstein copula weight matrix.
#' @param prior Numeric matrix (m × n). Dirichlet prior for Bernstein copula weight matrix.
#' @param options List. Control parameters. See Details.
#'
#' @return A list with elements:
#' \describe{
#'   \item{R}{Estimated weight matrix.}
#'   \item{llf}{Final log-likelihood.}
#'   \item{converged}{Logical. TRUE if convergence criteria met.}
#'   \item{iter}{Number of iterations.}
#'   \item{abserror}{Final absolute error.}
#'   \item{relerror}{Final relative error.}
#' }
#'
#' @details
#' The \code{options} list may include:
#' \itemize{
#'   \item \code{maxiter}: Maximum number of iterations (default: 1000)
#'   \item \code{abstol}: Absolute tolerance for log-likelihood (default: 1e-3)
#'   \item \code{reltol}: Relative tolerance (default: 1e-6)
#'   \item \code{mstep_maxiter}: Maximum M-step iterations (default: 1000)
#'   \item \code{mstep_abstol}: M-step tolerance (default: 1e-10)
#'   \item \code{mstep}: M-step algorithm ("sinkhorn" or "dou") (default: "sinkhorn")
#'   \item \code{verbose}: Logical. Print progress (default: FALSE)
#' }
#'
#' @export
#'
#' @examples
#' u <- runif(100)
#' v <- runif(100)
#' R0 <- matrix(1/9, 3, 3)
#' result <- emloop(u, v, R0)
#' result$R
emloop <- function(u, v, R_init, prior = NULL, options = list()) {
  stopifnot(is.matrix(R_init))
  m <- nrow(R_init)
  n <- ncol(R_init)
  stopifnot(length(u) == length(v))

  # default prior: ones
  if (is.null(prior)) {
    prior <- matrix(1.0, m, n)
  } else {
    stopifnot(all(dim(prior) == dim(R_init)))
  }

  defaults <- list(
    maxiter = 1000,
    abstol = 1e-3,
    reltol = 1e-6,
    mstep_maxiter = 1000,
    mstep_abstol = 1e-10,
    mstep = "sinkhorn", # or "dou"
    verbose = FALSE
  )
  opts <- modifyList(defaults, options)

  bernstein_emloop(u, v, R_init, prior, opts)
}
