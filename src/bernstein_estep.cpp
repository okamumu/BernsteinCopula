#include <RcppArmadillo.h>
using namespace Rcpp;

void dbinom_loader_inplace(int n, double p, NumericVector& v);

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double bernstein_estep(const arma::vec& u, const arma::vec& v,
  const arma::mat& R, const arma::mat& prior, arma::mat& tau) {
  int m = R.n_rows;
  int n = R.n_cols;
  int N = u.n_elem;

  NumericVector Bi(m);
  NumericVector Bj(n);

  arma::mat tmp(m, n);
  double llf = 0.0;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      tau(i, j) = prior(i, j);
    }
  }

  for (int k = 0; k < N; k++) {
    double uk = u(k);
    double vk = v(k);
    double total = 0.0;

    dbinom_loader_inplace(m-1, uk, Bi);
    dbinom_loader_inplace(n-1, vk, Bj);

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double val = R(i, j) * Bi(i) * Bj(j);
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
double bernstein_estep_with_weight(
  const arma::vec& u, const arma::vec& v,
  const arma::vec& du, const arma::vec& dv,
  const arma::mat& R,
  const arma::mat& prior,
  arma::mat& tau,
  arma::vec& exw1, arma::vec& exw2,
  arma::vec& eyw1, arma::vec& eyw2) {
  int m = R.n_rows;
  int n = R.n_cols;
  int N = u.n_elem;

  exw1.fill(0.0);
  exw2.fill(0.0);
  eyw1.fill(0.0);
  eyw2.fill(0.0);

  NumericVector Bi(m);
  NumericVector Bj(n);

  arma::mat tmp(m, n);
  double llf = 0.0;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      tau(i, j) = prior(i, j);
    }
  }

  for (int k = 0; k < N; k++) {
    double uk = u(k);
    double vk = v(k);
    double total = 0.0;

    dbinom_loader_inplace(m-1, uk, Bi);
    dbinom_loader_inplace(n-1, vk, Bj);

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double val = R(i, j) * Bi(i) * Bj(j);
        exw1(k) += i * val;
        exw2(k) += (m - i - 1)* val;
        eyw1(k) += j * val;
        eyw2(k) += (n - j - 1)* val;
        tmp(i, j) = val;
        total += val;
      }
    }

    if (total > 0) {
      llf += log(total) + log(du(k)) + log(dv(k));
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          tau(i, j) += tmp(i, j) / total;
        }
      }
      exw1(k) /= total;
      exw2(k) /= total;
      eyw1(k) /= total;
      eyw2(k) /= total;
    }
  }

  return llf;
}
