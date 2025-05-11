#include <RcppArmadillo.h>
using namespace Rcpp;

NumericVector dbinom_loader(int n, double p);

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double bernstein_copula_pdf(double u, double v,
                            const arma::mat& R) {
  int m = R.n_rows;
  int n = R.n_cols;

  NumericVector b_k = dbinom_loader(m-1, u);
  NumericVector b_l = dbinom_loader(n-1, v);

  double pdf = 0.0;
  for (int k = 0; k < m; ++k) {
    for (int l = 0; l < n; ++l) {
      double r_kl = R(k, l);
      pdf += r_kl * b_k(k) * b_l(l);
    }
  }

  return m * n * pdf;
}
