#include <RcppArmadillo.h>
using namespace Rcpp;

void dbinom_loader_inplace(int n, double p, NumericVector& v);

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double bernstein_estep(const arma::vec& u, const arma::vec& v, const arma::mat& R, arma::mat& tau) {
  int m = R.n_rows;
  int n = R.n_cols;
  int N = u.n_elem;

  NumericVector Bi(m);
  NumericVector Bj(n);

  arma::mat tmp(m, n);
  double llf = 0.0;

  for (int k = 0; k < N; k++) {
    double uk = u[k];
    double vk = v[k];
    double total = 0.0;

    dbinom_loader_inplace(m-1, uk, Bi);
    dbinom_loader_inplace(n-1, vk, Bj);

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double val = R(i, j) * Bi[i] * Bj[j];
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

// [[Rcpp::export]]
double bernstein_estep_with_weight(const arma::vec& u, const arma::vec& v, const arma::mat& R,
  arma::mat& tau, arma::vec& ew) {
  int m = R.n_rows;
  int n = R.n_cols;
  int N = u.n_elem;

  ew.fill(0.0);

  NumericVector Bi(m);
  NumericVector Bj(n);

  arma::mat tmp(m, n);
  double llf = 0.0;

  for (int k = 0; k < N; k++) {
    double uk = u[k];
    double vk = v[k];
    double total = 0.0;

    dbinom_loader_inplace(m-1, uk, Bi);
    dbinom_loader_inplace(n-1, vk, Bj);

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double val = R(i, j) * Bi[i] * Bj[j];
        ew[k] += i * val;
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
      ew[k] /= total;
    }
  }

  return llf;
}
