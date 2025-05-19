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
#' result <- fit.copula(u, v, R0)
#' result$R
fit.copula <- function(u, v, R_init, prior = NULL, options = list()) {
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
    maxiter = 2000,
    abstol = 1e-3,
    reltol = 1e-6,
    mstep_maxiter = 2000,
    mstep_abstol = 1e-10,
    mstep = "sinkhorn", # or "dou"
    verbose = FALSE
  )
  opts <- modifyList(defaults, options)

  bernstein_emloop(u, v, R_init, prior, opts)
}

#' EM Algorithm for Bernstein Copula with Flexible Marginals (via R6)
#'
#' Runs the EM algorithm to estimate a Bernstein copula and marginal distributions jointly.
#' Marginals are specified via R6 objects `Fx` and `Gy`, each of which should provide:
#' - `params`: named list of initial parameters,
#' - `pdf_func(x, params)`: function to evaluate the marginal density,
#' - `cdf_func(x, params)`: function to evaluate the marginal CDF,
#' - `emstep_func(x, w1, w2, params)`: function to update the parameters in the M-step.
#'
#' @param x Numeric vector. Samples for variable X.
#' @param y Numeric vector. Samples for variable Y.
#' @param R_init Matrix. Initial copula weight matrix of size \eqn{m \times n}.
#' @param prior Matrix or `NULL`. Prior weight matrix for tau-bar; defaults to a matrix of ones.
#' @param Fx R6 object representing the marginal model for X. Must have fields and methods as described above.
#' @param Gy R6 object representing the marginal model for Y. Must have fields and methods as described above.
#' @param options List of algorithm options. Must contain:
#'   \describe{
#'     \item{maxiter}{Maximum number of EM iterations.}
#'     \item{abstol}{Absolute tolerance for convergence.}
#'     \item{reltol}{Relative tolerance for convergence.}
#'     \item{joint.est}{Logical. If `TRUE`, update marginal parameters in the M-step.}
#'     \item{steps}{Number of steps for progress output.}
#'     \item{mstep_maxiter}{Maximum number of M-step iterations (for copula update).}
#'     \item{mstep_abstol}{Tolerance for M-step convergence.}
#'     \item{mstep}{Character string, either `"dou"` or `"sinkhorn"` to choose M-step method.}
#'     \item{verbose}{Logical. If `TRUE`, print progress during EM.}
#'   }
#'
#' @return A list containing:
#' \describe{
#'   \item{R}{Estimated copula weight matrix.}
#'   \item{params1}{Final parameters for the X marginal model.}
#'   \item{params2}{Final parameters for the Y marginal model.}
#'   \item{llf}{Final log-likelihood value.}
#'   \item{converged}{Logical. Whether the algorithm converged.}
#'   \item{iter}{Number of EM iterations completed.}
#'   \item{abserror}{Final absolute error in log-likelihood.}
#'   \item{relerror}{Final relative error in log-likelihood.}
#' }
#'
#' @export
#'
#' @examples
#' # See vignette or examples for full usage with R6 marginal models.
fit.copula.joint <- function(x, y, R_init, prior = NULL,
  Fx, Gy, options = list()) {
  stopifnot(is.matrix(R_init))
  m <- nrow(R_init)
  n <- ncol(R_init)
  stopifnot(length(x) == length(y))

  # default prior: ones
  if (is.null(prior)) {
    prior <- matrix(1.0, m, n)
  } else {
    stopifnot(all(dim(prior) == dim(R_init)))
  }

  defaults <- list(
    maxiter = 2000,
    abstol = 1e-3,
    reltol = 1e-6,
    joint.est = TRUE,
    steps = 100,
    mstep_maxiter = 2000,
    mstep_abstol = 1e-10,
    mstep = "sinkhorn", # or "dou"
    verbose = TRUE
  )
  opts <- modifyList(defaults, options)

  N <- length(x)
  exw1 <- numeric(N)
  exw2 <- numeric(N)
  eyw1 <- numeric(N)
  eyw2 <- numeric(N)

  R <- R_init
  params1 <- Fx$params
  params2 <- Gy$params
  llf_old <- -Inf

  converged <- FALSE
  for (iter in seq_len(opts$maxiter)) {

    tau_bar <- prior

    du <- Fx$pdf_func(x, params1)
    dv <- Gy$pdf_func(y, params2)
    u  <- Fx$cdf_func(x, params1)
    v  <- Gy$cdf_func(y, params2)

    llf <- bernstein_estep_with_weight(u, v, du, dv, R, tau_bar, exw1, exw2, eyw1, eyw2)

    if (opts$mstep == "dou") {
      R <- dou_mstep(tau_bar, maxiter = opts$mstep_maxiter, tol = opts$mstep_abstol)
    } else if (opts$mstep == "sinkhorn") {
      R <- sinkhorn_scaling(tau_bar, maxiter = opts$mstep_maxiter, tol = opts$mstep_abstol)
    } else {
      stop("Unknown mstep method: ", opts$mstep)
    }

    if (opts$joint.est == TRUE && !is.null(Fx$emstep_func) && !is.null(Gy$emstep_func)) {
      # Update marginal parameters
      params1 <- Fx$emstep_func(x, exw1, exw2, params1)
      params2 <- Gy$emstep_func(y, eyw1, eyw2, params2)
    }

    if (llf - llf_old < 0) {
      warning("Log-likelihood decreased: ", llf, " < ", llf_old)
    }
    abserror <- abs(llf - llf_old)
    relerror <- abs(llf - llf_old) / (abs(llf_old) + 1e-10)

    if (opts$verbose && iter %% opts$steps == 0) {
      cat(sprintf("iter: %d, loglik: %.6f, abs: %.2e, rel: %.2e\n", iter, llf, abserror, relerror))
    }

    if (abserror < opts$abstol && relerror < opts$reltol) {
      converged <- TRUE
      break
    }

    llf_old <- llf
  }

  list(
    R = R,
    params1 = params1,
    params2 = params2,
    llf = llf,
    converged = converged,
    iter = iter,
    abserror = abserror,
    relerror = relerror
  )
}