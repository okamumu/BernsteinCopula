#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double bisection(std::function<double(double)> f, double a, double b, double tol, int maxiter) {
  double fa = f(a), fb = f(b);
  if (fa * fb > 0) return NA_REAL;

  for (int i = 0; i < maxiter; ++i) {
    double c = 0.5 * (a + b);
    double fc = f(c);
    if (std::abs(fc) < tol || (b - a) < tol)
      return c;
    if (fa * fc < 0) {
      b = c; fb = fc;
    } else {
      a = c; fa = fc;
    }
  }
  return NA_REAL;
}

// [[Rcpp::export]]
arma::mat dou_mstep(const arma::mat& tau_bar, int maxiter = 1000, double tol = 1e-10) {
  int m = tau_bar.n_rows;
  int n = tau_bar.n_cols;

  arma::vec mu(m, arma::fill::value(0.5));
  arma::vec lambda(n, arma::fill::ones);
  arma::vec mu_new(m);

  for (int iter = 0; iter < maxiter; ++iter) {
    arma::vec lambda_old = lambda;
    arma::vec mu_old = mu;

    // M1: update lambda
    for (int l = 0; l < n; ++l) {
      double min_mu = mu.min();
      auto f_lambda = [&](double lam) {
        double sum = 0.0;
        for (int k = 0; k < m; ++k)
          sum += tau_bar(k, l) / (mu[k] + lam);
        return sum - 1.0 / n;
      };
      lambda[l] = bisection(f_lambda, -min_mu + 1e-8, 1e5, tol, 100);
    }

    // M2: update mu
    for (int k = 0; k < m; ++k) {
      double min_lambda = lambda.min();
      auto f_mu = [&](double muk) {
        double sum = 0.0;
        for (int l = 0; l < n; ++l)
          sum += tau_bar(k, l) / (muk + lambda[l]);
        return sum - 1.0 / m;
      };
      mu_new[k] = bisection(f_mu, -min_lambda + 1e-8, 1e5, tol, 100);
    }

    // M3: identifiability shift
    double shift = arma::mean(mu) - arma::mean(mu_new);
    mu = mu_new + shift;

    if (arma::norm(mu_new - mu, "inf") < tol &&
        arma::norm(lambda - lambda_old, "inf") < tol)
      break;
  }

  // Compute R
  arma::mat R(m, n);
  for (int k = 0; k < m; ++k)
    for (int l = 0; l < n; ++l)
      R(k, l) = tau_bar(k, l) / (mu[k] + lambda[l]);

  return R / arma::accu(R);
}
