% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em-group.R
\name{fit.copula.group.joint}
\alias{fit.copula.group.joint}
\title{EM Algorithm for Joint Bernstein Copula with Binned Data}
\usage{
fit.copula.group.joint(
  xbreaks,
  ybreaks,
  counts,
  R_init,
  prior = NULL,
  Fx,
  Gy,
  options = list()
)
}
\arguments{
\item{xbreaks}{Numeric vector. Breakpoints for the first variable.}

\item{ybreaks}{Numeric vector. Breakpoints for the second variable.}

\item{counts}{Numeric matrix. Counts of observations in each bin.}

\item{R_init}{Numeric matrix (m × n). Initial Bernstein copula weight matrix.}

\item{prior}{Numeric matrix (m × n). Dirichlet prior for Bernstein copula weight matrix.}

\item{Fx}{R6 object representing the marginal model for X. Must have fields and methods as described above.}

\item{Gy}{R6 object representing the marginal model for Y. Must have fields and methods as described above.}

\item{options}{List of algorithm options. Must contain:
\describe{
  \item{maxiter}{Maximum number of EM iterations.}
  \item{abstol}{Absolute tolerance for convergence.}
  \item{reltol}{Relative tolerance for convergence.}
  \item{joint.est}{Logical. If `TRUE`, update marginal parameters in the M-step.}
  \item{steps}{Number of steps for progress output.}
  \item{mstep_maxiter}{Maximum number of M-step iterations (for copula update).}
  \item{mstep_abstol}{Tolerance for M-step convergence.}
  \item{mstep}{Character string, either `"dou"` or `"sinkhorn"` to choose M-step method.}
  \item{verbose}{Logical. If `TRUE`, print progress during EM.}
}}
}
\value{
A list containing:
\describe{
  \item{R}{Estimated copula weight matrix.}
  \item{params1}{Final parameters for the X marginal model.}
  \item{params2}{Final parameters for the Y marginal model.}
  \item{llf}{Final log-likelihood value.}
  \item{converged}{Logical. Whether the algorithm converged.}
  \item{iter}{Number of EM iterations completed.}
  \item{abserror}{Final absolute error in log-likelihood.}
  \item{relerror}{Final relative error in log-likelihood.}
}
}
\description{
EM Algorithm for Joint Bernstein Copula with Binned Data
}
