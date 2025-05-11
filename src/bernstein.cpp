#include <RcppArmadillo.h>
using namespace Rcpp;

NumericVector dbinom_loader(int n, double p);

// [[Rcpp::depends(RcppArmadillo)]]

// Bernstein basis polynomial
inline double b_poly(int k, int n, double u) {
  return R::choose(n, k) * std::pow(u, k) * std::pow(1 - u, n - k);
}

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
      // double b_k = b_poly(k, m - 1, u);
      // double b_l = b_poly(l, n - 1, v);
      pdf += r_kl * b_k[k] * b_l[l];
    }
  }

  return m * n * pdf;
}
