#include <RcppArmadillo.h>
using namespace Rcpp;

void dbinom_loader_inplace(int n, double p, NumericVector& v);

// [[Rcpp::depends(RcppArmadillo)]]

// [0,1]x[0,1] area is divided into (n_x-1)x(n_y-1) rectangles where n_x = length(u), n_y = length(v)
// N is a (n_x-1)x(n_y-1) matrix meaning the number of samples observed in each area

// [[Rcpp::export]]
double bernstein_estep_group(const arma::vec& u, const arma::vec& v,
                             const arma::mat& R, const arma::mat& N,
                             const arma::mat& prior, arma::mat& tau) {
  int m = R.n_rows;
  int n = R.n_cols;
  int n_x = u.n_elem;
  int n_y = v.n_elem;

  arma::mat barF(m, n_x);
  arma::mat barG(n, n_y);
  NumericVector tmpF(m+1);
  NumericVector tmpG(n+1);

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      tau(i, j) = prior(i, j);
    }
  }

  // compute barF
  for (int k = 0; k < n_x; k++) {
    double uk = u(k);
    dbinom_loader_inplace(m, uk, tmpF);
    barF(0, k) = tmpF(0);
    for (int i = 1; i < m; i++) {
      barF(i, k) = barF(i-1, k) + tmpF(i);
    }
  }

  // compute barG
  for (int l = 0; l < n_y; l++) {
    double vl = v(l);
    dbinom_loader_inplace(n, vl, tmpG);
    barG(0, l) = tmpG(0);
    for (int j = 1; j < n; j++) {
      barG(j, l) = barG(j-1, l) + tmpG(j);
    }
  }

  // compute weights
  double llf = 0.0;
  for (int k = 1; k < n_x; k++) {
    for (int l = 1; l < n_y; l++) {
      double w = 0.0;
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          w += R(i, j) * (barF(i, k) - barF(i, k-1)) * (barG(j, l) - barG(j, l-1));
        }
      }
      if (w <= 0) {
        Rcpp::stop("Negative weight encountered in Bernstein E-step.");
      }
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          tau(i, j) += N(k-1, l-1) * R(i, j) * (barF(i, k) - barF(i, k-1)) * (barG(j, l) - barG(j, l-1)) / w;
        }
      }
      llf += N(k-1, l-1) * log(w); // - std::lgamma(N(k-1, l-1) + 1.0);
    }
  }

  return llf;
}

// [[Rcpp::export]]
double bernstein_estep_group_with_weight(
    const arma::vec& u, const arma::vec& v,
    const arma::mat& R,
    const arma::mat& N,
    const arma::mat& prior,
    arma::mat& tau,
    arma::vec& exw1, arma::vec& exw2,
    arma::vec& eyw1, arma::vec& eyw2) {
  int m = R.n_rows;
  int n = R.n_cols;
  int n_x = u.n_elem;
  int n_y = v.n_elem;

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      tau(i, j) = prior(i, j);
    }
  }
  exw1.fill(0.0);
  exw2.fill(0.0);
  eyw1.fill(0.0);
  eyw2.fill(0.0);

  arma::mat barF(m, n_x);
  arma::mat barFx1(m, n_x);
  arma::mat barFx2(m, n_x);
  arma::mat barG(n, n_y);
  arma::mat barGy1(n, n_y);
  arma::mat barGy2(n, n_y);
  NumericVector tmpF(m+1);
  NumericVector tmpG(n+1);

  // compute barF
  for (int k = 0; k < n_x; k++) {
    double uk = u(k);
    dbinom_loader_inplace(m, uk, tmpF);
    barF(0, k) = tmpF(0);
    barFx1(0, k) = 0.0;
    barFx2(0, k) = m * tmpF(0);
    for (int i = 1; i < m; i++) {
      barF(i, k) = barF(i-1, k) + tmpF(i);
      barFx1(i, k) = barFx1(i-1, k) + i * tmpF(i);
      barFx2(i, k) = barFx2(i-1, k) + (m - i) * tmpF(i);
    }
  }

  // compute barG
  for (int l = 0; l < n_y; l++) {
    double vl = v(l);
    dbinom_loader_inplace(n, vl, tmpG);
    barG(0, l) = tmpG(0);
    barGy1(0, l) = 0.0;
    barGy2(0, l) = n * tmpG(0);
    for (int j = 1; j < n; j++) {
      barG(j, l) = barG(j-1, l) + tmpG(j);
      barGy1(j, l) = barGy1(j-1, l) + j * tmpG(j);
      barGy2(j, l) = barGy2(j-1, l) + (n - j) * tmpG(j);
    }
  }

  // compute weights
  double llf = 0.0;
  for (int k = 1; k < n_x; k++) {
    for (int l = 1; l < n_y; l++) {
      if (N(k-1, l-1) == 0) {
        continue;
      }
      double w = 0.0;
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          w += R(i, j) * (barF(i, k) - barF(i, k-1)) * (barG(j, l) - barG(j, l-1));
        }
      }
      if (w <= 0) {
        Rcpp::stop("Negative weight encountered in Bernstein E-step. w=%e", w);
      }
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          tau(i, j) += N(k-1, l-1) * R(i, j) * (barF(i, k) - barF(i, k-1)) * (barG(j, l) - barG(j, l-1)) / w;
          exw1(k) += N(k-1, l-1) * R(i, j) * barFx1(i, k) * (barG(j, l) - barG(j, l-1)) / w;
          exw1(k-1) += -N(k-1, l-1) * R(i, j) * barFx1(i, k-1) * (barG(j, l) - barG(j, l-1)) / w;
          exw2(k) += N(k-1, l-1) * R(i, j) * barFx2(i, k) * (barG(j, l) - barG(j, l-1)) / w;
          exw2(k-1) += -N(k-1, l-1) * R(i, j) * barFx2(i, k-1) * (barG(j, l) - barG(j, l-1)) / w;
          eyw1(l) += N(k-1, l-1) * R(i, j) * (barF(i, k) - barF(i, k-1)) * barGy1(j, l) / w;
          eyw1(l-1) += -N(k-1, l-1) * R(i, j) * (barF(i, k) - barF(i, k-1)) * barGy1(j, l-1) / w;
          eyw2(l) += N(k-1, l-1) * R(i, j) * (barF(i, k) - barF(i, k-1)) * barGy2(j, l) / w;
          eyw2(l-1) += -N(k-1, l-1) * R(i, j) * (barF(i, k) - barF(i, k-1)) * barGy2(j, l-1) / w;
        }
      }
      llf += N(k-1, l-1) * log(w); // - std::lgamma(N(k-1, l-1) + 1.0);
    }
  }

  return llf;
}

