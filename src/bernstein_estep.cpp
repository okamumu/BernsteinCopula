#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double bernstein_em_estep(const arma::vec& u, const arma::vec& v, const arma::mat& R, arma::mat& tau) {
  int m = R.n_rows;
  int n = R.n_cols;
  int N = u.n_elem;

  arma::mat tmp(m, n);
  double llf = 0.0;

  for (int k = 0; k < N; k++) {
    double uk = u[k];
    double vk = v[k];
    double total = 0.0;

    for (int i = 0; i < m; i++) {
      double Bi = R::dbinom(i, m - 1, uk, false);
      for (int j = 0; j < n; j++) {
        double Bj = R::dbinom(j, n - 1, vk, false);
        double val = R(i, j) * Bi * Bj;
        tmp(i, j) = val;
        total += val;
      }
    }

    if (total > 0) {
      llf += log(total);
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          tau(i, j) += tmp(i, j) / total;
        }
      }
    }
  }

  return llf;
}
