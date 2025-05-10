#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat sinkhorn_scaling(const arma::mat& W, int maxiter = 1000, double tol = 1e-10) {
  int m = W.n_rows;
  int n = W.n_cols;

  arma::vec target_row(m, arma::fill::ones);
  arma::vec target_col(n, arma::fill::ones);
  target_row /= m;
  target_col /= n;

  arma::vec r(m, arma::fill::ones);
  arma::vec c(n, arma::fill::ones);

  for (int iter = 0; iter < maxiter; ++iter) {
    arma::vec r_old = r;
    arma::vec c_old = c;

    arma::vec Wc = W * c;
    for (int i = 0; i < m; ++i) {
      if (Wc[i] > 0) r[i] = target_row[i] / Wc[i];
    }

    arma::vec WTr = W.t() * r;
    for (int j = 0; j < n; ++j) {
      if (WTr[j] > 0) c[j] = target_col[j] / WTr[j];
    }

    if (arma::norm(r - r_old, "inf") < tol && arma::norm(c - c_old, "inf") < tol) {
      break;
    }
  }

  arma::mat R = arma::diagmat(r) * W * arma::diagmat(c);
  return R;
}
