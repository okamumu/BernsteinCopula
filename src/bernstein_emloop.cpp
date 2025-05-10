#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double bernstein_em_estep(const arma::vec& u, const arma::vec& v, const arma::mat& R, arma::mat& tau);
arma::mat dou_em_mstep(const arma::mat& tau_bar, int maxiter = 1000, double tol = 1e-10);
arma::mat sinkhorn_scaling(const arma::mat& W, int maxiter = 1000, double tol = 1e-10);

// [[Rcpp::export]]
List bernstein_emloop(const arma::vec& u,
  const arma::vec& v,
  arma::mat R_init,
  arma::mat tau_bar_init,
  List options) {
    // Extract options
    int maxiter = as<int>(options["maxiter"]);
    double abstol = as<double>(options["abstol"]);
    double reltol = as<double>(options["reltol"]);
    int mstep_maxiter = as<int>(options["mstep_maxiter"]);
    double mstep_abstol = as<double>(options["mstep_abstol"]);
    int verbose = as<bool>(options["verbose"]);
    std::string mstep_method = as<std::string>(options["mstep"]);

    int m = R_init.n_rows;
    int n = R_init.n_cols;

    arma::mat R = R_init;
    arma::mat tau_bar(m, n);
    double llf = 0.0, llf_old = -std::numeric_limits<double>::infinity();

    int iter = 0;
    bool converged = false;
    double abserror;
    double relerror;
    while (1) {
      tau_bar = tau_bar_init;
      llf = bernstein_em_estep(u, v, R, tau_bar);
      if (mstep_method == "dou") {
        R = dou_em_mstep(tau_bar, mstep_maxiter, mstep_abstol);
      } else if (mstep_method == "sinkhorn") {
        R = sinkhorn_scaling(tau_bar, mstep_maxiter, mstep_abstol);
      } else {
        stop("Unknown mstep method: " + mstep_method);
      }

      if (iter == 0) {
        llf_old = llf;
        ++iter;
        continue;
      }

      if (llf < llf_old) {
        warning("Log-likelihood decreased at iteration: %d (old: %g, new: %g)", iter, llf_old, llf);
      }

      abserror = std::abs(llf - llf_old);
      relerror = std::abs(llf - llf_old) / (std::abs(llf_old) + 1e-10);
      if (abserror < abstol && relerror < reltol) {
        converged = true;
        if (verbose) {
          Rcpp::Rcout << "Converged at iteration: " << iter << ", Log-likelihood: " << llf << std::endl;
        }
        break;
      }
      llf_old = llf;

      if (verbose) {
        Rcpp::Rcout << "Iteration: " << iter << ", Log-likelihood: " << llf << std::endl;
      }
      ++iter;
      Rcpp::checkUserInterrupt();

      if (iter >= maxiter) {
        if (verbose) {
          Rcpp::Rcout << "Maximum iterations reached: " << maxiter << std::endl;
        }
        break;
      }
    }

    return List::create(Named("R") = R,
    Named("llf") = llf,
    Named("converged") = converged,
    Named("iter") = iter,
    Named("abserror") = abserror,
    Named("relerror") = relerror);
  }
