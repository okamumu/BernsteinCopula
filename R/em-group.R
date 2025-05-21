#' EM Algorithm for Bernstein Copula with Binned Data
#'
#' @param ubreaks Numeric vector. Breakpoints for the first variable.
#' @param vbreaks Numeric vector. Breakpoints for the second variable.
#' @param counts Numeric matrix. Counts of observations in each bin.
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
#'   \item \code{maxiter}: Maximum number of iterations (default: 2000)
#'   \item \code{abstol}: Absolute tolerance for log-likelihood (default: 1e-3)
#'   \item \code{reltol}: Relative tolerance (default: 1e-6)
#'   \item \code{steps}: Number of steps for verbose output (default: 100)
#'   \item \code{mstep_maxiter}: Maximum M-step iterations (default: 2000)
#'   \item \code{mstep_abstol}: M-step tolerance (default: 1e-10)
#'   \item \code{mstep}: M-step algorithm ("sinkhorn" or "dou") (default: "sinkhorn")
#'   \item \code{verbose}: Logical. Print progress (default: FALSE)
#' }
#'
#' @export
#'
#' @examples
#' ubreaks <- seq(0, 1, length.out = 6)
#' vbreaks <- seq(0, 1, length.out = 6)
#' counts <- matrix(c(10, 20, 30, 40, 50,
#'                   15, 25, 35, 45, 55,
#'                   20, 30, 40, 50, 60,
#'                   25, 35, 45, 55, 65,
#'                   30, 40, 50, 60, 70), nrow = 5, byrow = TRUE)
#' R_init <- matrix(1/9, 3, 3)
#' result <- fit.copula.group(ubreaks, vbreaks, counts, R_init)
#' result$R
#' result$llf
#' result$converged
#' result$iter
#' result$abserror
#' result$relerror
fit.copula.group <- function(ubreaks, vbreaks, counts, R_init, prior = NULL, options = list()) {
  stopifnot(is.matrix(counts))
  stopifnot(is.matrix(R_init))
  m <- nrow(R_init)
  n <- ncol(R_init)

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
    steps = 100,
    mstep_maxiter = 2000,
    mstep_abstol = 1e-10,
    mstep = "sinkhorn", # or "dou"
    verbose = TRUE
  )
  opts <- modifyList(defaults, options)

  # run EM algorithm
  R <- R_init
  llf_old <- -Inf
  converged <- FALSE
  for (iter in seq_len(opts$maxiter)) {

    tau_bar <- matrix(0, m, n)
    llf <- bernstein_estep_group(ubreaks, vbreaks, R, counts, prior, tau_bar)

    if (opts$mstep == "dou") {
      R <- dou_mstep(tau_bar, maxiter = opts$mstep_maxiter, tol = opts$mstep_abstol)
    } else if (opts$mstep == "sinkhorn") {
      R <- sinkhorn_scaling(tau_bar, maxiter = opts$mstep_maxiter, tol = opts$mstep_abstol)
    } else if (opts$mstep == "noconstraint") {
      R <- tau_bar / sum(tau_bar)
    } else {
      stop("Unknown mstep method: ", opts$mstep)
    }

    if (llf - llf_old < 0) {
      warning("Log-likelihood decreased at iteration ", iter, ": ", llf, " < ", llf_old)
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
    llf = llf,
    converged = converged,
    iter = iter,
    abserror = abserror,
    relerror = relerror
  )
}

#' EM Algorithm for Joint Bernstein Copula with Binned Data
#'
#' @param xbreaks Numeric vector. Breakpoints for the first variable.
#' @param ybreaks Numeric vector. Breakpoints for the second variable.
#' @param counts Numeric matrix. Counts of observations in each bin.
#' @param R_init Numeric matrix (m × n). Initial Bernstein copula weight matrix.
#' @param prior Numeric matrix (m × n). Dirichlet prior for Bernstein copula weight matrix.
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
fit.copula.group.joint <- function(xbreaks, ybreaks, counts, R_init, prior = NULL,
  Fx, Gy, options = list()) {
  stopifnot(is.matrix(counts))
  stopifnot(is.matrix(R_init))
  m <- nrow(R_init)
  n <- ncol(R_init)

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
    steps = 100,
    joint.est = TRUE,
    mstep_maxiter = 2000,
    mstep_abstol = 1e-10,
    mstep = "sinkhorn", # or "dou"
    verbose = TRUE
  )
  opts <- modifyList(defaults, options)

  n_x <- length(xbreaks)
  n_y <- length(ybreaks)
  exw1 <- numeric(n_x)
  exw2 <- numeric(n_x)
  eyw1 <- numeric(n_y)
  eyw2 <- numeric(n_y)

  R <- R_init
  params1 <- Fx$params
  params2 <- Gy$params
  llf_old <- -Inf

  converged <- FALSE

  for (iter in seq_len(opts$maxiter)) {

    tau_bar <- matrix(0, m, n)

    u  <- Fx$cdf_func(xbreaks, params1)
    v  <- Gy$cdf_func(ybreaks, params2)

    llf <- bernstein_estep_group_with_weight(
      u, v, R, counts, prior, tau_bar, exw1, exw2, eyw1, eyw2
    )

    if (opts$mstep == "dou") {
      R <- dou_mstep(tau_bar, maxiter = opts$mstep_maxiter, tol = opts$mstep_abstol)
    } else if (opts$mstep == "sinkhorn") {
      R <- sinkhorn_scaling(tau_bar, maxiter = opts$mstep_maxiter, tol = opts$mstep_abstol)
    } else if (opts$mstep == "noconstraint") {
      R <- tau_bar / sum(tau_bar)
    } else {
      stop("Unknown mstep method: ", opts$mstep)
    }

    if (opts$joint.est == TRUE && !is.null(Fx$emstep_func) && !is.null(Gy$emstep_func)) {
      # Update marginal parameters
      params1 <- Fx$emstep_func_group(xbreaks, exw1, exw2, params1)
      params2 <- Gy$emstep_func_group(ybreaks, eyw1, eyw2, params2)
    }

    if (llf - llf_old < 0) {
      warning("Log-likelihood decreased at iteration ", iter, ": ", llf, " < ", llf_old)
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
